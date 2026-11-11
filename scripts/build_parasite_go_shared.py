'''
The host processes the parasites of one host reach, compared on equal terms. Every
parasite is tested against the host proteins both niches can reach, the extracellular
pool, on its interactors within that pool, so that a process can be tested for every
parasite whatever its niche. Species are averaged per genus, as species of one genus are
not independent, and a genus counts as enriched for a process at FDR <= 0.05 in at least
half its species. Groups are compared on the mean log2 odds ratios of their genera with a
two-sided Mann-Whitney test, Benjamini-Hochberg corrected over the processes enriched for
at least one genus.
Writes paper/tables/parasite_go_shared.csv (every contrast of every process) and
paper/figures/parasite_go_shared.svg/.png (processes by taxon).
'''
import argparse
import os
import sys
import textwrap

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [ROOT, os.path.join(ROOT, 'app')]

import utils
import web_utils
from scripts import figure_style
from scripts.build_parasite_go import GO_SEQUENTIAL, LABEL_CHARS, NICHE_COLORS

HELMINTHS = {'Nematoda', 'Trematoda', 'Cestoda'}
FLATWORMS = {'Trematoda', 'Cestoda'}
# rows of the figure per contrast direction
ROWS_PER_BLOCK = 6
# the share of genera of both helminths and unicellular parasites a process shared by
# the two has to be enriched for
SHARED_SHARE = 2 / 3
# log2 odds ratios are drawn up to this magnitude
COLOR_LIMIT = 4
WARM = ['#f6d9c2', '#e69f6b', '#c2601f', '#7a3408']


def log_odds(sets, background):
    '''log2 odds ratio of every testable term for every parasite, with a 0.5 correction.'''
    terms = utils.testable_terms(background)
    pairs = background['pairs']
    pairs = pairs[pairs['term'].isin(terms)]
    sizes = background['sizes'][terms]
    total = len(background['annotated'])
    odds = {}
    for taxid, proteins in sets.items():
        tested = proteins & background['annotated']
        a = (pairs[pairs['#string_protein_id'].isin(tested)].groupby('term').size()
             .reindex(terms, fill_value=0))
        b = len(tested) - a
        c = sizes - a
        d = total - sizes - b
        odds[taxid] = np.log2((a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5)))

    return pd.DataFrame(odds)


def genus_table(values, parasites):
    '''Parasite columns averaged per genus; a genus in two niches is two columns.'''
    units = {t: (parasites[t]['label'].split(' ')[0], parasites[t]['niche'],
                 parasites[t]['group']) for t in values.columns}
    grouped = values.T.groupby(pd.Series(units)).mean().T
    grouped.columns = pd.MultiIndex.from_tuples(grouped.columns,
                                                names=['genus', 'niche', 'group'])

    return grouped


def contrast(odds, tested_terms, first, second, name, labels):
    '''Mann-Whitney test of two sets of genera over the terms tested.'''
    rows = []
    for term in tested_terms:
        x, y = odds.loc[term, first], odds.loc[term, second]
        rows.append((name, labels[0], labels[1], term, x.mean(), y.mean(),
                     stats.mannwhitneyu(x, y, alternative='two-sided').pvalue))
    table = pd.DataFrame(rows, columns=['contrast', 'first', 'second', 'go_id', 'mean_first',
                                        'mean_second', 'p_value'])
    table['fdr_bh'] = multipletests(table['p_value'], method='fdr_bh')[1]

    return table


def draw(odds, enriched, blocks, names, columns, output_stem):
    figure_style.apply()
    labels, ys, spans = [], [], []
    y = 0.0
    wrapped = {term: textwrap.fill(names[term], LABEL_CHARS)
               for _, terms in blocks for term in terms}
    for title, terms in blocks:
        # one row of room for the title of the block
        y += 1
        positions = figure_style.row_positions([wrapped[t] for t in terms],
                                                figure_style.DENSE_ROW)
        spans.append((title, y - 1))
        ys += [y + p for p in positions]
        labels += [wrapped[t] for t in terms]
        y += positions[-1] + 1.5

    rows = [t for _, terms in blocks for t in terms]
    height = figure_style.DENSE_ROW * y + 2.15
    fig, ax = plt.subplots(figsize=(figure_style.WIDTH, height))
    cmap = LinearSegmentedColormap.from_list('diverging', WARM[::-1] + ['#f2f2f2'] + GO_SEQUENTIAL)
    norm = TwoSlopeNorm(vmin=-COLOR_LIMIT, vcenter=0, vmax=COLOR_LIMIT)

    xs = {}
    x = 0.0
    for i, column in enumerate(columns):
        if i and column['niche'] != columns[i - 1]['niche']:
            x += 0.25
        xs[i] = x
        x += 1

    for term, row_y in zip(rows, ys):
        for i, column in enumerate(columns):
            genera = column['genera']
            value = odds.loc[term, genera].mean()
            shown = np.clip(value, -COLOR_LIMIT, COLOR_LIMIT)
            ax.add_patch(Rectangle((xs[i] - 0.48, row_y - 0.46), 0.96, 0.92,
                                   color=cmap(norm(shown)), linewidth=0))
            hits = int(enriched.loc[term, genera].sum())
            if hits:
                ink = 'white' if abs(shown) > 0.55 * COLOR_LIMIT else '#222222'
                ax.text(xs[i], row_y, f'{hits}/{len(genera)}', ha='center', va='center',
                        color=ink, fontsize=figure_style.FONT)

    for title, row_y in spans:
        ax.text(-0.62, row_y, title, ha='right', va='center', fontweight='bold')

    ax.set_xlim(-0.6, x - 0.4)
    ax.set_ylim(y - 1.2, -0.6)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, linespacing=1.05)
    ax.set_xticks([xs[i] for i in range(len(columns))])
    ax.set_xticklabels([c['label'] for c in columns], rotation=90)
    ax.xaxis.tick_top()
    ax.tick_params(length=0, pad=2)
    for side in ax.spines.values():
        side.set_visible(False)

    # the niche over its columns, named on the strip
    strip = 1.3
    for niche, ink in (('intracellular', 'white'), ('extracellular', '#222222')):
        own = [xs[i] for i, c in enumerate(columns) if c['niche'] == niche]
        ax.add_patch(Rectangle((own[0] - 0.48, -0.7 - strip), own[-1] - own[0] + 0.96, strip,
                               color=NICHE_COLORS[niche], linewidth=0, clip_on=False))
        ax.text((own[0] + own[-1]) / 2, -0.7 - strip / 2, niche.capitalize(), ha='center',
                va='center', color=ink, clip_on=False)
    ax.tick_params(axis='x', pad=4 + (strip + 0.1) * figure_style.DENSE_ROW * 72)

    fig.subplots_adjust(left=0.40, right=0.985, top=1 - 1.45 / height, bottom=0.62 / height)
    key = fig.add_axes([0.40, 0.3 / height, 0.3, 0.08 / height])
    points = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    colorbar = fig.colorbar(points, cax=key, orientation='horizontal', extend='both')
    colorbar.set_label('Mean log$_2$ odds ratio over the genera', labelpad=2)
    colorbar.ax.xaxis.set_label_position('top')
    colorbar.outline.set_visible(False)
    for extension in ('svg', 'png'):
        fig.savefig(f'{output_stem}.{extension}', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--host', type=int, default=9606, help='host taxid (default: human)')
    parser.add_argument('--min-score', type=float, default=0.35,
                        help='links below this confidence score are left out')
    parser.add_argument('--fdr', type=float, default=0.05)
    parser.add_argument('--table', default='paper/tables/parasite_go_shared.csv')
    parser.add_argument('--figure', default='paper/figures/parasite_go_shared',
                        help='output path without extension')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    parasites = config['parasites']
    predictions = pd.read_parquet(os.path.join(args.data_dir, 'predictions.parquet'))
    edges = predictions[(predictions['taxid2'] == str(args.host))
                        & (predictions['weight'].astype(float) >= args.min_score)]
    go_df = web_utils.go_background(args.data_dir, (str(args.host),), 'extracellular')
    background = utils.prepare_go_background(go_df)
    # a parasite with no interactor in the pool has nothing to compare
    sets = {int(t): set(rows['target']) & background['annotated']
            for t, rows in edges.groupby('taxid1')}
    sets = {t: s for t, s in sets.items() if s}

    names = background['names']
    odds = genus_table(log_odds(sets, background), parasites)
    tests = pd.concat(utils.calculate_enrichment(s, go_df, background=background).assign(taxid=t)
                      for t, s in sets.items())
    hits = tests[tests['fdr_bh'] <= args.fdr].assign(hit=True)
    enriched = (hits.pivot_table(index='go_id', columns='taxid', values='hit', aggfunc='any')
                .reindex(index=odds.index, columns=list(sets)).notna())
    genus_hits = genus_table(enriched.astype(float), parasites) >= 0.5
    tested_terms = genus_hits.index[genus_hits.any(axis=1)]

    def genera(test):
        return [c for c in odds.columns if test(c)]

    helminths = genera(lambda c: c[2] in HELMINTHS)
    unicellular = genera(lambda c: c[2] not in HELMINTHS)
    flatworms = genera(lambda c: c[2] in FLATWORMS)
    nematodes = genera(lambda c: c[2] == 'Nematoda' and c[1] == 'extracellular')
    table = pd.concat([
        contrast(odds, tested_terms, flatworms, nematodes,
                 'flatworms vs extracellular nematodes', ('flatworms', 'nematodes')),
        contrast(odds, tested_terms, helminths, unicellular,
                 'helminths vs unicellular parasites', ('helminths', 'unicellular')),
        contrast(odds, tested_terms, genera(lambda c: c[1] == 'extracellular'),
                 genera(lambda c: c[1] == 'intracellular'),
                 'extracellular vs intracellular', ('extracellular', 'intracellular')),
        contrast(odds, tested_terms,
                 genera(lambda c: c[1] == 'extracellular' and c[2] not in HELMINTHS),
                 genera(lambda c: c[1] == 'intracellular' and c[2] not in HELMINTHS),
                 'unicellular: extracellular vs intracellular', ('extracellular', 'intracellular')),
    ], ignore_index=True)
    table['go_term'] = table['go_id'].map(names)
    table['difference'] = table['mean_first'] - table['mean_second']
    summary = table[table['fdr_bh'] <= args.fdr].groupby('contrast')['difference']
    print(pd.DataFrame({'higher in first': summary.apply(lambda d: (d > 0).sum()),
                        'higher in second': summary.apply(lambda d: (d < 0).sum())}))
    table.sort_values(['contrast', 'p_value']).to_csv(args.table, index=False)
    print(f'Wrote {args.table}: {len(tested_terms)} processes in {table["contrast"].nunique()} contrasts')

    # the rows: most significant first, within each direction, leaving out a process whose
    # host proteins largely repeat one already shown
    members = utils.go_term_members([go_df])
    picked = []

    def pick(candidates, limit=ROWS_PER_BLOCK):
        chosen = []
        for term in candidates:
            if utils.representing_term(term, picked, members) is None:
                picked.append(term)
                chosen.append(term)
            if len(chosen) == limit:
                break
        return chosen

    def direction(name, sign, higher):
        '''Terms of a contrast higher on one side and enriched for half its genera.'''
        rows = table[(table['contrast'] == name) & (table['fdr_bh'] <= args.fdr)
                     & (np.sign(table['difference']) == sign)]
        common = genus_hits[higher].mean(axis=1) >= 0.5
        rows = rows[rows['go_id'].map(common)]
        return rows.sort_values(['p_value', 'go_id'])['go_id'].tolist()

    shares = pd.DataFrame({'helminths': genus_hits[helminths].mean(axis=1),
                           'unicellular': genus_hits[unicellular].mean(axis=1)})
    shared = shares[shares.min(axis=1) >= SHARED_SHARE]
    blocks = [
        ('Flatworms over nematodes',
         pick(direction('flatworms vs extracellular nematodes', 1, flatworms))),
        ('Nematodes over flatworms',
         pick(direction('flatworms vs extracellular nematodes', -1, nematodes))),
        ('Helminths over unicellular',
         pick(direction('helminths vs unicellular parasites', 1, helminths))),
        ('Unicellular over helminths',
         pick(direction('helminths vs unicellular parasites', -1, unicellular))),
        ('Shared by both', pick(shared.min(axis=1).sort_values(ascending=False).index.tolist())),
    ]
    blocks = [block for block in blocks if block[1]]

    columns = []
    for niche in NICHE_COLORS:
        for group in config['parasite_groups']:
            own = genera(lambda c: c[1] == niche and c[2] == group)
            if own:
                columns.append({'niche': niche, 'genera': own,
                                'label': f'{group} ({len(own)})'})
    draw(odds, genus_hits, blocks, names, columns, args.figure)
    print(f'Wrote {args.figure}.svg and .png: {sum(len(b[1]) for b in blocks)} processes')
