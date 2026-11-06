'''
Which host orthologous groups the parasites of each niche converge on, ranked within the
niche so that the 26 extracellular parasites do not outvote the 17 intracellular ones.
Groups are counted in genera, as species of one genus are not independent: a genus reaches
a group when at least half its species of the niche do. A group counts by the genera that
still reach it without its largest parasite group, the one giving it the most links, so
that one conserved parasite family cannot place a group on its own.
Writes paper/tables/shared_ogs_niche.csv (per niche, every group reached by two or more of
its parasites, with its largest parasite group and how many genera of the other niche
reach it) and paper/figures/shared_ogs_niche.svg/.png (the groups that --min-intracellular
and --min-extracellular genera reach without their largest parasite group, as one dot
matrix per niche over that niche's parasites).
'''
import argparse
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils
from scripts import figure_style
from scripts.build_shared_ogs import DOT_COLOR, row_label

NICHES = ['intracellular', 'extracellular']
# inches: the labels left of the columns, one column, the header over the columns, the
# gap between the panels, and the keys at the foot
GUTTER = 2.0
COLUMN = 0.155
HEADER = 1.2
GAP = 0.2
LEGEND = 0.75
# columns the margin counts start after, and take each
MARGIN_START = 0.6
MARGIN_STEP = 1.3


def genera_reached(edges, species):
    '''Per (niche, group), the genera of which at least half the species of the niche reach the group.'''
    reached = edges.groupby(['niche', 'group2', 'genus'])['taxid1'].nunique().reset_index()
    needed = pd.MultiIndex.from_frame(reached[['niche', 'genus']]).map(species)
    reached = reached[reached['taxid1'] * 2 >= needed]

    return reached.groupby(['niche', 'group2']).size()


def count_niche_ogs(config, predictions, host):
    '''One row per (niche, host orthologous group) reached by two or more parasites of the niche.'''
    parasites = config['parasites']
    edges = predictions[predictions['taxid2'] == str(host)].copy()
    edges['niche'] = edges['taxid1'].map(lambda t: parasites[int(t)]['niche'])
    edges['genus'] = edges['taxid1'].map(lambda t: parasites[int(t)]['label'].split()[0])
    species = edges.groupby(['niche', 'genus'])['taxid1'].nunique().to_dict()

    # the parasite group giving each group the most links
    links = edges.groupby(['niche', 'group2', 'group1']).size().rename('links').reset_index()
    links = links.sort_values(['links', 'group1'], ascending=[False, True], kind='stable')
    largest = links.drop_duplicates(['niche', 'group2']).set_index(['niche', 'group2'])['group1']
    edge_keys = pd.MultiIndex.from_frame(edges[['niche', 'group2']])
    rest = edges[largest.reindex(edge_keys).to_numpy() != edges['group1'].to_numpy()]

    genera = genera_reached(edges, species)
    left = genera_reached(rest, species)
    ogs = edges.groupby(['niche', 'group2']).agg(
        parasites=('taxid1', 'nunique'),
        parasite_groups=('group1', 'nunique'),
        host_proteins=('target', 'nunique'),
        interactions=('target', 'size'),
        symbols=('target_name', lambda n: ', '.join(sorted(set(str(x).upper() for x in n))))).reset_index()
    other = {'intracellular': 'extracellular', 'extracellular': 'intracellular'}
    keys = list(zip(ogs['niche'], ogs['group2']))
    ogs.insert(3, 'genera', [genera.get(k, 0) for k in keys])
    ogs.insert(4, 'genera_without_largest', [left.get(k, 0) for k in keys])
    ogs.insert(5, 'largest_parasite_group', [largest[k] for k in keys])
    ogs.insert(6, 'other_niche_genera', [genera.get((other[n], g), 0) for n, g in keys])
    ogs = ogs[ogs['parasites'] > 1].rename(columns={'group2': 'orthologous_group'})
    ogs['niche'] = pd.Categorical(ogs['niche'], NICHES)
    ogs = ogs.sort_values(['niche', 'genera_without_largest', 'genera', 'parasites', 'parasite_groups', 'interactions'],
                          ascending=[True, False, False, False, False, False], kind='stable')

    return ogs, edges


def draw_panel(fig, config, ogs, edges, niche, minimum, top, size):
    '''The matrix of one niche, its top `top` inches below the top of the figure. Returns its bottom.'''
    parasites = config['parasites']
    taxon_order = list(config['parasite_groups'])
    width, height = fig.get_size_inches()
    columns = sorted(set(edges.loc[edges['niche'] == niche, 'taxid1']),
                     key=lambda t: (taxon_order.index(parasites[int(t)]['group']),
                                    parasites[int(t)]['label']))
    col_of = {t: i for i, t in enumerate(columns)}
    labels = [row_label(g, s) for g, s in zip(ogs['orthologous_group'], ogs['symbols'])]
    ys = figure_style.row_positions(labels)
    row_of = dict(zip(ogs['orthologous_group'], ys))
    n_cols = len(columns)
    n_genera = edges.loc[edges['niche'] == niche, 'genus'].nunique()

    ax_width = (n_cols + 0.2) * COLUMN
    ax_height = (ys[-1] + 1.2) * figure_style.ROW
    bottom = top + HEADER + ax_height
    ax = fig.add_axes([GUTTER / width, 1 - bottom / height, ax_width / width, ax_height / height])

    # a dot per (group, parasite), its area the parasite proteins reaching the group
    dots = (edges[(edges['niche'] == niche) & edges['group2'].isin(row_of)]
            .groupby(['group2', 'taxid1'])['source'].nunique().reset_index())
    ax.scatter([col_of[t] for t in dots['taxid1']], [row_of[g] for g in dots['group2']],
               s=dots['source'].map(size), color=DOT_COLOR, linewidths=0, zorder=3)
    ax.set_xlim(-0.6, n_cols - 0.4)
    ax.set_ylim(ys[-1] + 0.6, -0.6)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, linespacing=1.05)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([figure_style.short_name(parasites[int(t)]['label']) for t in columns],
                       rotation=90, style='italic')
    ax.tick_params(length=0, pad=2)
    ax.grid(color='#eeeeee', linewidth=0.5)
    ax.set_axisbelow(True)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.xaxis.tick_top()

    # the taxonomic strip over the columns, under the parasite names
    strip = 0.55
    for i, t in enumerate(columns):
        ax.add_patch(Rectangle((i - 0.5, -0.6 - strip), 1, strip * 0.8, linewidth=0, clip_on=False,
                               color=config['parasite_groups'][parasites[int(t)]['group']]))
    ax.tick_params(axis='x', pad=2 + strip * figure_style.ROW * 72)

    headers = [('genera', 'genera'), ('genera_without_largest', 'without largest'),
               ('other_niche_genera', 'other niche'),
               ('parasite_groups', 'parasite groups')]
    for j, (column, header) in enumerate(headers):
        x = n_cols - 0.4 + MARGIN_START + MARGIN_STEP * j
        ax.text(x, -0.6 - strip - 0.3, header, ha='center', va='bottom', rotation=90,
                color='#555555', clip_on=False)
        for row in ogs.itertuples():
            ax.text(x, row_of[row.orthologous_group], f'{getattr(row, column)}', ha='center',
                    va='center', color='#555555', clip_on=False)

    fig.text(0.06 / width, 1 - (top + HEADER - 0.02) / height,
             f'{niche.capitalize()} parasites\nreached by ≥{minimum} of {n_genera} genera\nwithout the largest parasite group',
             ha='left', va='bottom', fontsize=figure_style.FONT, linespacing=1.15, color='#333333')

    return bottom


def draw(config, panels, edges, output_stem):
    figure_style.apply()
    size = lambda n: 6 + 3 * min(n, 20)
    heights = [HEADER + (figure_style.row_positions(
        [row_label(g, s) for g, s in zip(ogs['orthologous_group'], ogs['symbols'])])[-1] + 1.2)
        * figure_style.ROW for _, ogs, _ in panels]
    fig = plt.figure(figsize=(figure_style.WIDTH, sum(heights) + GAP + LEGEND))

    top = 0.05
    for niche, ogs, minimum in panels:
        top = draw_panel(fig, config, ogs, edges, niche, minimum, top, size) + GAP

    # the size key sits beside the taxa, as in the GO figure
    groups = [Patch(color=color, label=group) for group, color in config['parasite_groups'].items()
              if group in {config['parasites'][int(t)]['group'] for t in set(edges['taxid1'])}]
    sizes = [Line2D([], [], linestyle='', marker='o', color=DOT_COLOR, markersize=size(n) ** 0.5,
                    label=label) for n, label in ((1, '1'), (5, '5'), (20, '20+'))]
    figure_style.stack_legends(fig, left=0.01, blocks=[('Taxonomic group', groups)])
    fig.legend(handles=sizes, title='Parasite proteins', loc='lower left', ncol=len(sizes),
               frameon=False, handlelength=0.9, columnspacing=1.2, alignment='left',
               bbox_to_anchor=(0.705, 0.02 / fig.get_figheight()), title_fontsize=figure_style.FONT)
    for extension in ('svg', 'png'):
        fig.savefig(f'{output_stem}.{extension}', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--host', type=int, default=9606, help='host taxid (default: human)')
    parser.add_argument('--min-intracellular', type=int, default=5,
                        help='groups fewer intracellular genera reach without their largest parasite group stay out of that panel')
    parser.add_argument('--min-extracellular', type=int, default=13,
                        help='groups fewer extracellular genera reach without their largest parasite group stay out of that panel')
    parser.add_argument('--table', default='paper/tables/shared_ogs_niche.csv')
    parser.add_argument('--figure', default='paper/figures/shared_ogs_niche',
                        help='output path without extension')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    predictions = pd.read_parquet(os.path.join(args.data_dir, 'predictions.parquet'))
    ogs, edges = count_niche_ogs(config, predictions, args.host)
    ogs.to_csv(args.table, index=False)
    minimums = {'intracellular': args.min_intracellular, 'extracellular': args.min_extracellular}
    panels = []
    for niche in NICHES:
        kept = ogs[(ogs['niche'] == niche) & (ogs['genera_without_largest'] >= minimums[niche])]
        panels.append((niche, kept, minimums[niche]))
        print(f"{niche}: {(ogs['niche'] == niche).sum()} groups reached by two or more parasites, "
              f"{len(kept)} by {minimums[niche]} or more genera without the largest parasite group, "
              f"{(kept['other_niche_genera'] == 0).sum()} of them by no genus of the other niche")
    print(f'Wrote {args.table}')
    draw(config, panels, edges, args.figure)
    print(f'Wrote {args.figure}.svg and .png')
