'''
The host processes the parasites of one host reach. For each parasite, the GO biological
processes over-represented among its predicted host interactors are tested as the app
tests them: one-sided Fisher's exact test against the host proteins of the parasite's
niche the pipeline passed (web_utils.go_background), Benjamini-Hochberg corrected. Each
parasite contributes its most significant processes as rows, skipping a process whose
host proteins largely repeat one already shown, and the rows are clustered by the
parasites they are enriched for.
Writes paper/tables/parasite_go.csv (every process enriched for every parasite, and
whether it is a row of the figure) and paper/figures/parasite_go.svg/.png (the dot
matrix).
'''
import argparse
import os
import sys
import textwrap

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [ROOT, os.path.join(ROOT, 'app')]

import utils
import web_utils
from scripts import figure_style

NICHE_COLORS = {'intracellular': '#3d3d3d', 'extracellular': '#c7c7c7'}
# the colour ramp of the app's enrichment figures
GO_SEQUENTIAL = ['#bcdcec', '#7fc0dd', '#3f9fca', '#2b8cbe', '#12587d', '#08324a']
# two processes are one row when this share of their host proteins is common to both
REDUNDANT_JACCARD = 0.5
# process names wrap onto lines of at most this many characters
LABEL_CHARS = 40
# the left margin, a figure fraction, that the wrapped names take
LEFT = 0.37


def enrich(config, data_dir, predictions, host, min_score):
    '''Every process tested for every parasite of the host, one row per parasite and term.'''
    edges = predictions[(predictions['taxid2'] == str(host))
                        & (predictions['weight'].astype(float) >= min_score)]
    backgrounds, tests = {}, []
    for taxid, rows in edges.groupby('taxid1'):
        niche = web_utils.parasite_niche(config, taxid)
        if niche not in backgrounds:
            backgrounds[niche] = web_utils.go_background(data_dir, (str(host),), niche)
        test = utils.calculate_enrichment(set(rows['target']), backgrounds[niche])
        if not test.empty:
            tests.append(test.drop(columns=['B', 'C', 'D'])
                         .rename(columns={'A': 'n_proteins'})
                         .assign(taxid=int(taxid), niche=niche))

    return pd.concat(tests, ignore_index=True), backgrounds


def pick_rows(enriched, backgrounds, per_parasite):
    '''
    Each parasite's most significant processes, a process standing in for any later one
    whose host proteins in the background it largely shares.
    '''
    members = {}
    for go_df in backgrounds.values():
        for term, ids in go_df.groupby('term')['#string_protein_id']:
            members.setdefault(term, set()).update(ids)

    def jaccard(a, b):
        return len(members[a] & members[b]) / len(members[a] | members[b])

    rows = []
    for _, terms in enriched.sort_values(['taxid', 'fdr_bh', 'go_id']).groupby('taxid'):
        covered = set()
        for term in terms['go_id']:
            match = next((r for r in rows if r == term
                          or jaccard(r, term) >= REDUNDANT_JACCARD), None)
            if match is None:
                rows.append(term)
                match = term
            covered.add(match)
            if len(covered) >= per_parasite:
                break

    return rows


def order_rows(enriched, rows, columns):
    '''The rows clustered by the parasites they are enriched for.'''
    presence = (enriched[enriched['go_id'].isin(rows)].assign(hit=1)
                .pivot_table(index='go_id', columns='taxid', values='hit', fill_value=0)
                .reindex(index=rows, columns=columns, fill_value=0))
    if len(rows) < 3:
        return rows
    distances = pdist(presence.to_numpy().astype(bool), metric='jaccard')
    tree = hierarchy.optimal_leaf_ordering(hierarchy.linkage(distances, 'average'), distances)

    return [rows[i] for i in hierarchy.leaves_list(tree)]


def row_positions(labels):
    '''
    The y of each row, in rows of figure_style.ROW: a row is one ROW from the next
    unless their wrapped names need more room.
    '''
    line = 1.05 * figure_style.FONT / 72 / figure_style.ROW
    halves = [line * (label.count('\n') + 1) / 2 for label in labels]
    ys = [0.0]
    for above, below in zip(halves, halves[1:]):
        ys.append(ys[-1] + max(1, above + below + 0.1))

    return ys


def draw(config, enriched, rows, columns, output_stem):
    figure_style.apply()
    parasites = config['parasites']
    taxon_colors = config['parasite_groups']
    col_of = {t: i for i, t in enumerate(columns)}
    names = dict(zip(enriched['go_id'], enriched['go_term']))

    dots = enriched[enriched['go_id'].isin(rows)].copy()
    odds = dots['odds_ratio'].astype(float)
    cap = odds[np.isfinite(odds)].max()
    # an infinite ratio is drawn at the largest finite one
    dots['log_odds'] = np.log2(odds.where(np.isfinite(odds), cap).clip(lower=1))

    labels = [textwrap.fill(names[r], LABEL_CHARS) for r in rows]
    ys = row_positions(labels)
    y_of = dict(zip(rows, ys))
    n_cols = len(columns)
    fig, ax = plt.subplots(figsize=(figure_style.WIDTH, figure_style.ROW * (ys[-1] + 1) + 2.35))
    cmap = LinearSegmentedColormap.from_list('go', GO_SEQUENTIAL)
    norm = Normalize(vmin=0, vmax=np.ceil(dots['log_odds'].max()))
    # the area of a dot grows with the proteins annotated
    size = lambda n: 4 + 1.1 * n
    points = ax.scatter([col_of[t] for t in dots['taxid']], [y_of[r] for r in dots['go_id']],
                        s=size(dots['n_proteins'].clip(upper=60)), c=dots['log_odds'],
                        cmap=cmap, norm=norm, linewidths=0, zorder=3)
    ax.set_xlim(-0.6, n_cols - 0.4)
    ax.set_ylim(ys[-1] + 0.6, -0.6)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, linespacing=1.05)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([figure_style.short_name(parasites[t]['label']) for t in columns],
                       rotation=90, style='italic')
    ax.tick_params(length=0, pad=2)
    ax.grid(color='#eeeeee', linewidth=0.5)
    ax.set_axisbelow(True)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.xaxis.tick_top()

    # two strips over the columns: niche, then taxonomic group
    strip = 0.55
    for i, t in enumerate(columns):
        p = parasites[t]
        ax.add_patch(Rectangle((i - 0.5, -0.6 - 2 * strip), 1, strip * 0.8,
                               color=NICHE_COLORS[p['niche']], linewidth=0, clip_on=False))
        ax.add_patch(Rectangle((i - 0.5, -0.6 - strip), 1, strip * 0.8,
                               color=taxon_colors[p['group']], linewidth=0, clip_on=False))
    ax.tick_params(axis='x', pad=2 + 2 * strip * figure_style.ROW * 72)

    groups = [Patch(color=color, label=group) for group, color in taxon_colors.items()
              if any(parasites[t]['group'] == group for t in columns)]
    niches = [Patch(color=color, label=niche) for niche, color in NICHE_COLORS.items()]
    sizes = [Line2D([], [], linestyle='', marker='o', color='#7fc0dd',
                    markersize=np.sqrt(size(n)), label=str(n)) for n in (5, 20, 60)]
    legends = figure_style.stack_legends(fig, left=LEFT, blocks=[
        ('Taxonomic group', groups), ('Niche', niches),
        ('Host proteins annotated to the process', sizes)])
    height = fig.get_figheight()
    fig.subplots_adjust(left=LEFT, right=0.98, top=1 - 1.2 / height,
                        bottom=(legends + 0.1) / height)
    # the colour key sits beside the size key, at the foot of the figure
    key = fig.add_axes([0.7, 0.2 / height, 0.2, 0.07 / height])
    colorbar = fig.colorbar(points, cax=key, orientation='horizontal')
    colorbar.set_label('log$_2$ odds ratio', labelpad=2)
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
    parser.add_argument('--per-parasite', type=int, default=2,
                        help='distinct processes each parasite contributes as rows')
    parser.add_argument('--table', default='paper/tables/parasite_go.csv')
    parser.add_argument('--figure', default='paper/figures/parasite_go',
                        help='output path without extension')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    predictions = pd.read_parquet(os.path.join(args.data_dir, 'predictions.parquet'))
    tested, backgrounds = enrich(config, args.data_dir, predictions, args.host,
                                 args.min_score)
    enriched = tested[tested['fdr_bh'] <= args.fdr]
    parasites = config['parasites']
    # by niche, then taxonomic group, then name, as the similarity figure orders them
    columns = sorted(enriched['taxid'].unique(),
                     key=lambda t: (list(NICHE_COLORS).index(parasites[t]['niche']),
                                    list(config['parasite_groups']).index(parasites[t]['group']),
                                    parasites[t]['label']))
    rows = order_rows(enriched, pick_rows(enriched, backgrounds, args.per_parasite), columns)

    table = enriched.assign(parasite=enriched['taxid'].map(lambda t: parasites[t]['label']),
                            in_figure=enriched['go_id'].isin(rows))
    table[['parasite', 'taxid', 'niche', 'go_id', 'go_term', 'n_proteins', 'odds_ratio',
           'p_value', 'fdr_bh', 'in_figure', 'nodes']].to_csv(args.table, index=False)
    print(f"Wrote {args.table}: {len(table)} enriched processes over {len(columns)} "
          f"parasites, {table['go_id'].nunique()} distinct")
    draw(config, enriched, rows, columns, args.figure)
    print(f"Wrote {args.figure}.svg and .png: {len(rows)} processes")
