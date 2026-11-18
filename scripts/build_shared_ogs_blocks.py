'''
The host orthologous groups the niches share and those each reaches mostly on its own, in
one heatmap over the parasites of one host. The shared block holds the groups reached by at
least --min-share of the parasites of each niche (build_shared_ogs). The niche blocks hold
the other groups reached by at least --min-extracellular or --min-intracellular genera
without their largest parasite group (build_shared_ogs_niche). A niche whose host pool holds
no protein of a group cannot reach it, and is hatched.
Writes paper/tables/shared_ogs_blocks.csv (the rows of the figure) and
paper/figures/shared_ogs_blocks.svg/.png (the heatmap).
'''
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils
from scripts.build_shared_ogs import LOC_CLASSES, LOC_NONE, LOC_SEVERAL
from scripts import build_shared_ogs, figure_style
from scripts.build_shared_ogs_niche import count_niche_ogs


# the DeepLoc 2 Accurate-model cut-offs the pipeline filters with (app/web_utils.py)
CUTOFFS = {'Extracellular': 0.61728516, 'Cell membrane': 0.56464844,
           'Cytoplasm': 0.47612305, 'Nucleus': 0.50136719}
DEEPLOC_COLUMNS = dict(zip(LOC_CLASSES, ('extracellular', 'cell_membrane', 'cytoplasm', 'nucleus')))


def unreachable_niches(data_dir, host, groups):
    '''{group: the niches whose host pool holds none of its proteins in this host}.'''
    eligible = pd.read_parquet(os.path.join(data_dir, 'eligible_proteins.parquet'))
    eligible = eligible[eligible['taxid'].astype(str) == str(host)].set_index('protein')
    orthologs = pd.read_parquet(os.path.join(data_dir, 'host_orthologs.parquet'))
    orthologs = orthologs[(orthologs['taxid'].astype(str) == str(host)) & orthologs['group'].isin(groups)]
    members = {row.group: row.proteins.split(',') for row in orthologs.itertuples()}

    return {group: {niche for niche in ('intracellular', 'extracellular')
                    if not eligible.reindex(members.get(group, []))[f'reachable_{niche}'].eq(True).any()}
            for group in groups}


def reached_localisation(data_dir, edges, groups):
    '''
    {group: class}: the DeepLoc class most of the reached host proteins of the group are
    called for, as the app reads it. A protein called for more than one class is Several,
    and so is a group on a tie.
    '''
    deeploc = pd.read_parquet(os.path.join(data_dir, 'deeploc_localisations.parquet')).set_index('protein')
    called = pd.DataFrame({c: deeploc[column] > CUTOFFS[c] for c, column in DEEPLOC_COLUMNS.items()})
    per_protein = called.apply(lambda row: LOC_SEVERAL if row.sum() > 1
                               else row.idxmax() if row.any() else LOC_NONE, axis=1)
    out = {}
    for group, targets in edges[edges['group2'].isin(groups)].groupby('group2')['target']:
        tally = per_protein.reindex(sorted(set(targets))).dropna().value_counts()
        if tally.empty:
            continue
        out[group] = (LOC_SEVERAL if len(tally) > 1 and tally.iloc[0] == tally.iloc[1]
                      else tally.index[0])

    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--host', type=int, default=9606, help='host taxid (default: human)')
    parser.add_argument('--min-share', type=float, default=0.5,
                        help='the share of the parasites of each niche that makes a group shared')
    parser.add_argument('--min-extracellular', type=int, default=13,
                        help='genera that must reach a group without its largest parasite group')
    parser.add_argument('--min-intracellular', type=int, default=5,
                        help='genera that must reach a group without its largest parasite group')
    parser.add_argument('--table', default='paper/tables/shared_ogs_blocks.csv')
    parser.add_argument('--figure', default='paper/figures/shared_ogs_blocks',
                        help='output path without extension')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    predictions = pd.read_parquet(os.path.join(args.data_dir, 'predictions.parquet'))
    ogs, edges = build_shared_ogs.count_ogs(config, predictions, args.host)
    by_niche, _ = count_niche_ogs(config, predictions, args.host)

    shared = list(build_shared_ogs.shared_ogs(ogs, edges, args.min_share)['orthologous_group'])
    blocks = [('Shared by the niches', shared)]
    for niche, minimum in (('extracellular', args.min_extracellular), ('intracellular', args.min_intracellular)):
        own = by_niche[(by_niche['niche'] == niche) & (by_niche['genera_without_largest'] >= minimum)]
        blocks.append((f'Mainly {niche}', [g for g in own['orthologous_group'] if g not in shared]))
    rows = [group for _, groups in blocks for group in groups]
    unreachable = unreachable_niches(args.data_dir, args.host, rows)

    drawn = ogs.set_index('orthologous_group').loc[rows].reset_index()
    genera = by_niche.pivot(index='orthologous_group', columns='niche', values='genera_without_largest')
    table = drawn.assign(
        block=[name for name, groups in blocks for _ in groups],
        family=drawn['orthologous_group'].map(lambda g: build_shared_ogs.GROUP_NAMES.get(g, '')),
        extracellular_genera_without_largest=genera['extracellular'].reindex(rows).fillna(0).astype(int).to_numpy(),
        intracellular_genera_without_largest=genera['intracellular'].reindex(rows).fillna(0).astype(int).to_numpy(),
        unreachable=[', '.join(sorted(unreachable[g])) for g in rows])
    table[['block', 'orthologous_group', 'family', 'parasites', 'intracellular', 'extracellular',
           'intracellular_genera_without_largest', 'extracellular_genera_without_largest',
           'parasite_groups', 'unreachable', 'symbols']].to_csv(args.table, index=False)
    print(f'Wrote {args.table}: ' + ', '.join(f'{len(g)} {name}' for name, g in blocks)
          + f", {sum('extracellular' in u for u in unreachable.values())} unreachable by extracellular parasites")
    build_shared_ogs.draw(config, drawn, edges, args.figure, row=figure_style.DENSE_ROW,
                          blocks=[(name, len(groups)) for name, groups in blocks], unreachable=unreachable,
                          localisation=reached_localisation(args.data_dir, edges, rows))
    print(f'Wrote {args.figure}.svg and .png: {len(rows)} groups')
