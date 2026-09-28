'''
How many STRING proteins of each host TISSUES annotates to each tissue in the config, at
one score cutoff for every host, as a supplementary figure. Writes
paper/tables/tissues_coverage.csv and paper/figures/tissues_coverage.svg/.png.
'''
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import utils
from scripts import figure_style
import compare_tissue_filters
import venn_tissues_vs_hpa as venn

BAR_COLOR = '#2a62ab'
CUTOFF = 1.0


def count(config_file, cutoff):
    '''One row per host and tissue: the host's STRING proteins annotated to it.'''
    hosts = utils.read_config(filepath=config_file, field='hosts')
    pools = compare_tissue_filters.host_pools(config_file)
    rows = []
    for taxid, host in hosts.items():
        by_tissue = venn.raw_tissues_tissue_proteins(config_file, pools[taxid], taxid, cutoff)
        for n, tissue in sorted(((len(p), t) for t, p in by_tissue.items()), reverse=True):
            rows.append({'host': host['label'], 'proteome': len(pools[taxid]),
                         'tissue': tissue, 'proteins': n})

    return pd.DataFrame(rows)


def draw(table, cutoff, output_stem):
    figure_style.apply()
    hosts = list(dict.fromkeys(table['host']))
    sizes = table.groupby('host').size()
    # hosts two to a row, each row as tall as its longer panel at the shared row height
    pairs = [hosts[i:i + 2] for i in range(0, len(hosts), 2)]
    heights = [max(sizes[h] for h in pair) for pair in pairs]
    head, gap, foot = 0.45, 0.45, 0.4
    fig_height = figure_style.ROW * sum(heights) + head * len(pairs) + gap * (len(pairs) - 1) + foot
    fig = plt.figure(figsize=(figure_style.WIDTH, fig_height))
    grid = fig.add_gridspec(len(pairs), 2, height_ratios=heights, wspace=0.55,
                            hspace=(head + gap) / (figure_style.ROW * sum(heights) / len(pairs)),
                            top=1 - head / fig_height, bottom=foot / fig_height, left=0.14, right=0.97)
    top = table['proteins'].max()
    for r, pair in enumerate(pairs):
        for c, host in enumerate(pair):
            ax = fig.add_subplot(grid[r, c])
            rows = table[table['host'] == host].reset_index(drop=True)
            ax.barh(range(len(rows)), rows['proteins'], color=BAR_COLOR, height=0.72)
            for i, n in enumerate(rows['proteins']):
                ax.text(n + 0.01 * top, i, f'{n:,}', va='center', color='#555555')
            ax.set_yticks(range(len(rows)))
            ax.set_yticklabels(rows['tissue'])
            # the shared row height, so bars are as thick in every panel
            ax.set_ylim(heights[r] - 0.4, -0.6)
            ax.set_xlim(0, top * 1.15)
            ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(5000))
            ax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))
            ax.tick_params(length=0)
            ax.grid(axis='x', color='#e6e6e6', linewidth=0.6)
            ax.set_axisbelow(True)
            for side in ('top', 'right', 'left'):
                ax.spines[side].set_visible(False)
            ax.spines['bottom'].set_color('#999999')
            letter = 'ABCDEFGH'[2 * r + c]
            ax.set_title(f'{letter}  {host}', loc='left', fontweight='bold', pad=14)
            ax.annotate(f"{len(rows)} tissues, {rows['proteome'].iloc[0]:,} STRING proteins",
                        (0, 1), xycoords='axes fraction', xytext=(0, 4), textcoords='offset points',
                        color='#555555')
            if r == len(pairs) - 1:
                ax.set_xlabel(f'Proteins with TISSUES score ≥ {cutoff:g}')
    for extension in ('svg', 'png'):
        fig.savefig(f'{output_stem}.{extension}', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--cutoff', type=float, default=CUTOFF, help='TISSUES score cutoff for every host')
    parser.add_argument('--table', default='paper/tables/tissues_coverage.csv')
    parser.add_argument('--figure', default='paper/figures/tissues_coverage', help='output path without extension')
    args = parser.parse_args()

    table = count(args.config, args.cutoff)
    table.to_csv(args.table, index=False)
    print(f'{len(table)} host tissues; wrote {args.table}')
    draw(table, args.cutoff, args.figure)
    print(f'Wrote {args.figure}.svg and .png')
