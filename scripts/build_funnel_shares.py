'''
Draw the share of proteins each stage of the pipeline keeps from the stage before it,
per host-parasite pair, from the counts build_filter_funnel.py writes. The funnel
figure's log axis shows fold reductions, so a late stage with a small count can look
wide; here every cell is kept / previous on the same 0-100% scale. One row per pair,
in the funnel's order and blocks, the host stages left and the parasite stages right.
Writes paper/tables/funnel_shares.csv and paper/figures/funnel_shares.svg/.png.
'''
import argparse
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils
from scripts import figure_style
from scripts.build_filter_funnel import NICHE_COLORS, NICHE_ORDER

# (column header, stage kept, stage before it), host side then parasite side
HOST_STEPS = [('TISSUES', 'host_tissue', 'host_proteome'),
              ('DeepLoc', 'host_localization', 'host_tissue'),
              ('KOG', 'host_kog', 'host_localization'),
              ('PPI', 'host_interactor', 'host_kog')]
PARASITE_STEPS = [('DeepLoc', 'parasite_localization', 'parasite_proteome'),
                  ('KOG', 'parasite_kog', 'parasite_localization'),
                  ('PPI', 'parasite_interactor', 'parasite_kog')]
# the scale of the similarity figure, white at nothing kept
SCALE = ['#ffffff', '#deebf7', '#9ecae1', '#6baed6', '#3182bd', '#08519c']
EMPTY_COLOR = '#d9d9d9'
# cells darker than this share take white text
DARK = 0.6
# the layout in inches: a cell is square-ish and wide enough for "100"; the parasite
# columns sit SIDE_GAP cells right of the host ones; the strips sit left of the cells
CELL = 0.27
SIDE_GAP = 0.45
STRIP_X = (-0.8, -0.5)
STRIP_WIDTH = 0.22
NAMES = 1.15
LEFT, TOP, BOTTOM, COLUMN_GAP = 0.05, 1.05, 0.3, 0.15
# rows between the stacked blocks of the right column, for the next header
BLOCK_GAP = 7.5


def shares(table):
    '''The share of each step per pair, empty where the stage before kept nothing.'''
    out = table[['host', 'species', 'group', 'niche']].copy()
    for side, steps in (('host', HOST_STEPS), ('parasite', PARASITE_STEPS)):
        for label, kept, before in steps:
            out[f'{side} {label}'] = (table[kept] / table[before].replace(0, np.nan)).round(4)

    return out


def x_positions():
    '''The left edge of every cell column in cell units, host then parasite.'''
    host = list(range(len(HOST_STEPS)))
    parasite = [len(HOST_STEPS) + SIDE_GAP + i for i in range(len(PARASITE_STEPS))]

    return host, parasite


def draw_block(ax, rows, group_colors, cmap, host_name):
    '''One host's pairs: a row of cells per pair, the share printed in each, the
    niche and group strips left of the cells and the parasite names left of those.'''
    host_x, parasite_x = x_positions()
    columns = [(x, f'host {label}') for x, (label, _, _) in zip(host_x, HOST_STEPS)] + \
              [(x, f'parasite {label}') for x, (label, _, _) in zip(parasite_x, PARASITE_STEPS)]
    for i, row in rows.iterrows():
        for x, column in columns:
            value = row[column]
            empty = pd.isna(value)
            ax.add_patch(Rectangle((x, i), 1, 1, facecolor=EMPTY_COLOR if empty else cmap(value),
                                   edgecolor='white', linewidth=1))
            if not empty:
                text = '<1' if 0 < value < 0.005 else f'{value * 100:.0f}'
                ax.text(x + 0.5, i + 0.52, text, ha='center', va='center',
                        color='white' if value > DARK else '#333333')
        # niche outermost, taxonomic group beside it, as in the funnel
        for x, color in zip(STRIP_X, (NICHE_COLORS[row['niche']], group_colors[row['group']])):
            ax.add_patch(Rectangle((x, i + 0.04), STRIP_WIDTH, 0.92, color=color, linewidth=0))

    right = parasite_x[-1] + 1
    ax.set_xlim(STRIP_X[0] - 0.05, right)
    ax.set_ylim(len(rows), 0)
    ax.set_yticks([i + 0.5 for i in range(len(rows))])
    ax.set_yticklabels([figure_style.short_name(s) for s in rows['species']], fontstyle='italic')
    ax.tick_params(length=0, pad=2)
    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # the step names over the cells, the side over them, the host over all
    for x, (label, _, _) in [*zip(host_x, HOST_STEPS), *zip(parasite_x, PARASITE_STEPS)]:
        ax.text(x + 0.5, -0.3, label, ha='center', va='bottom', rotation=90)
    for xs, side in ((host_x, 'Host'), (parasite_x, 'Parasite')):
        ax.annotate('', xy=(xs[0] + 0.05, -4.3), xytext=(xs[-1] + 0.95, -4.3),
                    arrowprops=dict(arrowstyle='-', color='#999999', linewidth=0.6),
                    annotation_clip=False)
        ax.text((xs[0] + xs[-1] + 1) / 2, -4.5, side, ha='center', va='bottom')
    ax.text(STRIP_X[0] - 0.05, -4.5, host_name, ha='right', va='bottom', fontweight='bold',
            fontsize=figure_style.TITLE)


def draw(table, group_colors, output_stem):
    '''The host with the most pairs in the left column, the other hosts stacked in the
    right, with the colour bar and legends under them.'''
    figure_style.apply()
    cmap = LinearSegmentedColormap.from_list('kept', SCALE)
    hosts = list(dict.fromkeys(table['host']))
    sizes = {h: (table['host'] == h).sum() for h in hosts}
    first = max(hosts, key=sizes.get)
    others = [h for h in hosts if h != first]
    row = figure_style.ROW
    height = sizes[first] * row + TOP + BOTTOM
    fig = plt.figure(figsize=(figure_style.WIDTH, height))
    _, parasite_x = x_positions()
    width = (parasite_x[-1] + 1 - STRIP_X[0] + 0.05) * CELL
    figure_width = figure_style.WIDTH

    def add_axes(x, y_top, n):
        return fig.add_axes([x / figure_width, (y_top - n * row) / height,
                             width / figure_width, n * row / height])

    x_first = LEFT + NAMES
    draw_block(add_axes(x_first, height - TOP, sizes[first]),
               table[table['host'] == first].reset_index(drop=True), group_colors, cmap, first)
    x_other = x_first + width + COLUMN_GAP + NAMES
    y = height - TOP
    for host in others:
        draw_block(add_axes(x_other, y, sizes[host]),
                   table[table['host'] == host].reset_index(drop=True), group_colors, cmap, host)
        y -= (sizes[host] + BLOCK_GAP) * row

    # the colour bar under the right column, the legends under it
    bar = fig.add_axes([(x_other - NAMES + 0.1) / figure_width, (y + BLOCK_GAP * row - 0.9) / height,
                        (width + NAMES - 0.4) / figure_width, 0.1 / height])
    colorbar = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap), cax=bar, orientation='horizontal')
    colorbar.set_ticks([0, 0.25, 0.5, 0.75, 1], labels=['0', '25', '50', '75', '100%'])
    colorbar.outline.set_visible(False)
    bar.set_title('Proteins kept from the stage before (%)', fontsize=figure_style.FONT, loc='left')
    figure_style.stack_legends(fig, left=(x_other - NAMES) / figure_width, ncol=3, blocks=[
        ('Parasite niche in the host',
         [Patch(color=NICHE_COLORS[niche], label=niche) for niche in NICHE_ORDER]),
        ('Taxonomic group', [Patch(color=color, label=group) for group, color in group_colors.items()])])
    for extension in ('svg', 'png'):
        fig.savefig(f'{output_stem}.{extension}', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--counts', default='paper/tables/funnel.csv')
    parser.add_argument('--table', default='paper/tables/funnel_shares.csv')
    parser.add_argument('--figure', default='paper/figures/funnel_shares', help='output path without extension')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    table = shares(pd.read_csv(args.counts))
    table.to_csv(args.table, index=False)
    print(f"Wrote {args.table} ({len(table)} pairs)")
    draw(table, config['parasite_groups'], args.figure)
    print(f"Wrote {args.figure} (.svg and .png)")
