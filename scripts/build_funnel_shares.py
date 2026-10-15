'''
Draw the share of proteins each stage of the pipeline keeps from the stage before it,
per host-parasite pair, from the counts build_filter_funnel.py writes. Every share is
on the same 0-100% scale, so a filter that removes much reads as a pale cell whatever
the size of the proteome. Two counts frame the shares: the parasite's STRING proteome
before its stages and the predicted interactions after them; the host's proteome is one
number per host, given in the block title. One row per pair, in the funnel's order and
blocks. Writes paper/tables/funnel_shares.csv and paper/figures/funnel_shares.svg/.png.
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
COUNT_COLOR = '#555555'
# cells darker than this share take white text
DARK = 0.6
# the layout in inches; a cell is at most CELL wide, narrower when the two columns of
# blocks would not otherwise fit the page. In cell units: the gap between the host and
# parasite columns, the one before the interactions, and the strips left of the cells
CELL = 0.27
SIDE_GAP = 0.45
COUNT_GAP = 0.3
STRIP_X = (-0.8, -0.5)
STRIP_WIDTH = 0.22
NAMES = 1.0
LEFT, RIGHT, TOP, BOTTOM, COLUMN_GAP = 0.05, 0.12, 1.25, 0.3, 0.15
# rows between the stacked blocks of the right column, for the next header
BLOCK_GAP = 9


def shares(table):
    '''Per pair: the share each stage keeps of the stage before, empty where that kept
    nothing, and the counts that frame them.'''
    out = table[['host', 'species', 'group', 'niche', 'host_proteome', 'parasite_proteome']].copy()
    for side, steps in (('host', HOST_STEPS), ('parasite', PARASITE_STEPS)):
        for label, kept, before in steps:
            out[f'{side} {label}'] = (table[kept] / table[before].replace(0, np.nan)).round(4)
    out['interactions'] = table['edges']

    return out


def layout():
    '''[(left edge in cell units, column, header, is a count)] of one row of cells,
    with the edges of the host and parasite groups for the labels over them.'''
    cells = [(i, f'host {label}', label, False) for i, (label, _, _) in enumerate(HOST_STEPS)]
    host = (0, len(HOST_STEPS))
    start = len(HOST_STEPS) + SIDE_GAP
    cells.append((start, 'parasite_proteome', 'Proteome', True))
    cells += [(start + 1 + i, f'parasite {label}', label, False)
              for i, (label, _, _) in enumerate(PARASITE_STEPS)]
    parasite = (start, start + 1 + len(PARASITE_STEPS))
    cells.append((parasite[1] + COUNT_GAP, 'interactions', 'Interactions', True))

    return cells, {'Host': host, 'Parasite': parasite}


def count_label(value):
    '''A count in at most four characters: thousands abbreviated.'''
    if value >= 9950:
        return f'{value / 1000:.0f}k'

    return f'{value / 1000:.1f}k' if value >= 1000 else f'{value:.0f}'


def draw_block(ax, rows, group_colors, cmap):
    '''One host's pairs: a row of cells per pair, the share printed in each and the
    counts in plain text, the niche and group strips left of the cells and the parasite
    names left of those.'''
    cells, groups = layout()
    for i, row in rows.iterrows():
        for x, column, _, is_count in cells:
            value = row[column]
            if is_count:
                ax.text(x + 0.5, i + 0.52, count_label(value), ha='center', va='center',
                        color=COUNT_COLOR)
                continue
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

    ax.set_xlim(STRIP_X[0] - 0.05, cells[-1][0] + 1)
    ax.set_ylim(len(rows), 0)
    ax.set_yticks([i + 0.5 for i in range(len(rows))])
    ax.set_yticklabels([figure_style.short_name(s) for s in rows['species']], fontstyle='italic')
    ax.tick_params(length=0, pad=2)
    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # the column names over the cells, the side over its group
    for x, _, header, is_count in cells:
        ax.text(x + 0.5, -0.3, header, ha='center', va='bottom', rotation=90,
                color=COUNT_COLOR if is_count else 'black')
    for side, (start, end) in groups.items():
        ax.annotate('', xy=(start + 0.05, -4.3), xytext=(end - 0.05, -4.3),
                    arrowprops=dict(arrowstyle='-', color='#999999', linewidth=0.6),
                    annotation_clip=False)
        ax.text((start + end) / 2, -4.5, side, ha='center', va='bottom')


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
    cells, _ = layout()
    units = cells[-1][0] + 1 - STRIP_X[0] + 0.05
    cell = min(CELL, (figure_style.WIDTH - LEFT - RIGHT - 2 * NAMES - COLUMN_GAP) / 2 / units)
    width = units * cell
    figure_width = figure_style.WIDTH

    def block(host, x, y_top):
        rows = table[table['host'] == host].reset_index(drop=True)
        ax = fig.add_axes([x / figure_width, (y_top - len(rows) * row) / height,
                           width / figure_width, len(rows) * row / height])
        draw_block(ax, rows, group_colors, cmap)
        # the host and its proteome on a line of their own over the side labels
        fig.text((x - NAMES) / figure_width, (y_top + 6.2 * row) / height,
                 f"{host} ({rows['host_proteome'].iloc[0]:,} proteins)", ha='left', va='bottom',
                 fontweight='bold', fontsize=figure_style.TITLE)

    x_first = LEFT + NAMES
    block(first, x_first, height - TOP)
    x_other = x_first + width + COLUMN_GAP + NAMES
    y = height - TOP
    for host in others:
        block(host, x_other, y)
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
