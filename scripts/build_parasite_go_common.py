'''
The host processes most parasites reach, selected as the shared groups of the human
interactome are. Every process is tested for every parasite as build_parasite_go does, and
the rows are chosen by how many parasites are enriched for them, in three blocks: the
processes enriched in at least --min-share of the parasites of each niche, then the
commonest others of each niche, --per-niche each. A process annotated to more than
--max-term host proteins stays out, as common and uninformative, and one that largely
repeats a row already chosen is left out as it is in build_parasite_go.
Writes paper/tables/parasite_go_common.csv (the rows, with the parasites and genera of each
niche enriched for them), paper/tables/parasite_go_enriched.csv (every process enriched for
every parasite) and paper/figures/parasite_go_common.svg/.png (the dot matrix).
'''
import argparse
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [ROOT, os.path.join(ROOT, 'app')]

import utils
from scripts import build_parasite_go
from scripts.build_parasite_go import NICHE_COLORS

NICHES = list(NICHE_COLORS)


def count_terms(enriched, members, config, taxids):
    '''One row per process: the parasites and genera of each niche enriched for it, its size and median odds.'''
    parasites = config['parasites']
    enriched = enriched.assign(genus=enriched['taxid'].map(lambda t: parasites[t]['label'].split()[0]))
    # the species of each genus among `taxids`, the parasites of this host
    species = {}
    for taxid in taxids:
        key = (parasites[taxid]['niche'], parasites[taxid]['label'].split()[0])
        species[key] = species.get(key, 0) + 1
    counts = enriched.groupby(['go_id', 'niche'])['taxid'].nunique().unstack(fill_value=0)
    counts = counts.reindex(columns=NICHES, fill_value=0)
    # a genus is enriched when at least half its species are
    by_genus = enriched.groupby(['go_id', 'niche', 'genus'])['taxid'].nunique().reset_index()
    needed = pd.MultiIndex.from_frame(by_genus[['niche', 'genus']]).map(species)
    genera = (by_genus[by_genus['taxid'] * 2 >= needed].groupby(['go_id', 'niche']).size()
              .unstack(fill_value=0).reindex(index=counts.index, columns=NICHES, fill_value=0))
    odds = enriched['odds_ratio'].astype(float).replace(np.inf, np.nan)

    return counts.assign(
        intracellular_genera=genera['intracellular'],
        extracellular_genera=genera['extracellular'],
        size=counts.index.map(lambda term: len(members[term])),
        median_odds=odds.groupby(enriched['go_id']).median(),
        go_term=enriched.drop_duplicates('go_id').set_index('go_id')['go_term'])


def pick_rows(counts, sizes, members, share, per_niche, max_term):
    '''The rows of each block, commonest first and the higher odds ratio on a tie: [(block, [term])].'''
    counts = counts[counts['size'] <= max_term]
    common = {niche: counts[niche] >= share * sizes[niche] for niche in NICHES}
    both = common['intracellular'] & common['extracellular']
    total = counts['intracellular'] + counts['extracellular']
    order = counts.assign(total=total).sort_values(['total', 'median_odds'], ascending=False)

    picked = utils.first_distinct_terms(list(order.index[both.reindex(order.index)]), members, len(order))
    blocks = [('Shared by the niches', list(picked))]
    for niche in reversed(NICHES):
        block = []
        for term in counts[common[niche] & ~both].sort_values([niche, 'median_odds'], ascending=False).index:
            if len(block) == per_niche:
                break
            if utils.representing_term(term, picked, members) is None:
                picked.append(term)
                block.append(term)
        blocks.append((f'Mainly {niche}', block))

    return blocks


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--host', type=int, default=9606, help='host taxid (default: human)')
    parser.add_argument('--min-score', type=float, default=0.35,
                        help='links below this confidence score are left out')
    parser.add_argument('--fdr', type=float, default=0.05)
    parser.add_argument('--min-share', type=float, default=0.5,
                        help='the share of the parasites of a niche that makes a process common in it')
    parser.add_argument('--per-niche', type=int, default=12,
                        help='rows of the commonest processes of each niche outside the shared ones')
    parser.add_argument('--max-term', type=int, default=200,
                        help='processes annotated to more host proteins stay out of the figure')
    parser.add_argument('--table', default='paper/tables/parasite_go_common.csv')
    parser.add_argument('--enriched-table', default='paper/tables/parasite_go_enriched.csv',
                        help='every process enriched for every parasite')
    parser.add_argument('--figure', default='paper/figures/parasite_go_common',
                        help='output path without extension')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    parasites = config['parasites']
    predictions = pd.read_parquet(os.path.join(args.data_dir, 'predictions.parquet'))
    tested, backgrounds = build_parasite_go.enrich(config, args.data_dir, predictions, args.host,
                                                   args.min_score)
    enriched = tested[tested['fdr_bh'] <= args.fdr]
    members = utils.go_term_members(backgrounds.values())
    # every parasite of the host is a column, as in the shared groups figure
    taxids = {int(t) for t in predictions.loc[predictions['taxid2'] == str(args.host), 'taxid1']}
    sizes = pd.Series([parasites[t]['niche'] for t in taxids]).value_counts()
    columns = sorted(taxids, key=lambda t: (NICHES.index(parasites[t]['niche']),
                                            list(config['parasite_groups']).index(parasites[t]['group']),
                                            parasites[t]['label']))

    counts = count_terms(enriched, members, config, taxids)
    blocks = pick_rows(counts, sizes, members, args.min_share, args.per_niche, args.max_term)
    rows = [term for _, terms in blocks for term in terms]
    table = pd.concat([counts.loc[terms].assign(block=name) for name, terms in blocks])
    table.reset_index()[['block', 'go_id', 'go_term', 'intracellular', 'extracellular',
                         'intracellular_genera', 'extracellular_genera', 'size',
                         'median_odds']].to_csv(args.table, index=False)
    print(f'Wrote {args.table}: ' + ', '.join(f'{len(terms)} {name}' for name, terms in blocks))
    (enriched.assign(parasite=enriched['taxid'].map(lambda t: parasites[t]['label']))
     [['parasite', 'taxid', 'niche', 'go_id', 'go_term', 'n_proteins', 'odds_ratio',
       'p_value', 'fdr_bh', 'nodes']].to_csv(args.enriched_table, index=False))
    build_parasite_go.draw(config, enriched, rows, columns, backgrounds, args.figure,
                           blocks=[(name, len(terms)) for name, terms in blocks])
    print(f'Wrote {args.figure}.svg and .png: {len(rows)} processes')
