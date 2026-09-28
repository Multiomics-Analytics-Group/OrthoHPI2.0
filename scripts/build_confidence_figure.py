'''
The confidence scores of the predicted interactions: one block per host, as tall as the
number of parasites infecting it, a box per parasite over the scores of its interactions,
laid out as the networks figure so the two are read together. A score is the average of
the experimental and database evidence of the KOG-KOG link the prediction was transferred
from, and a link is kept only if one of the two is at least 0.7, so no score is under
0.35. Writes paper/tables/confidence.csv (one row per host-parasite pair) and
paper/figures/confidence.svg/.png.
'''
import argparse
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils
from pipeline import homology
from scripts import figure_style

NICHE_COLORS = {'extracellular': '#c7c7c7', 'intracellular': '#3d3d3d'}
# a host with fewer parasites still gets a block this many boxes tall
MIN_COLUMN = 3
# the lowest score a kept link can have
SCORE_FLOOR = homology.EVIDENCE_CUTOFF / 2


def collect(config, data_dir):
    '''The score of every prediction, labelled with its host and parasite; hosts in config
    order, parasites by group, niche then name.'''
    parasites = config['parasites']
    hosts = {str(t): h['label'].split('(')[1].rstrip(')') for t, h in config['hosts'].items()}
    predictions = pd.read_parquet(os.path.join(data_dir, 'predictions.parquet'))
    predictions['weight'] = predictions['weight'].astype(float)
    order = list(config['parasite_groups'])

    frames = []
    for host, host_name in hosts.items():
        for taxid, parasite in sorted(parasites.items(),
                                      key=lambda item: (order.index(item[1]['group']),
                                                        list(NICHE_COLORS).index(item[1]['niche']),
                                                        item[1]['label'])):
            if int(host) not in parasite['hosts']:
                continue
            edges = predictions[(predictions['taxid1'] == str(taxid)) & (predictions['taxid2'] == host)]
            frames.append(pd.DataFrame({'host': host_name, 'parasite': parasite['label'],
                                        'group': parasite['group'], 'niche': parasite['niche'],
                                        'score': edges['weight'].values}))

    return pd.concat(frames, ignore_index=True)


def summarise(scores):
    '''One row per host-parasite pair: the interactions and the quartiles of their scores.'''
    pairs = scores.groupby(['host', 'parasite', 'group', 'niche'], sort=False)['score']
    table = pairs.describe()[['count', 'min', '25%', '50%', '75%', 'max']]
    table.columns = ['interactions', 'min', 'q1', 'median', 'q3', 'max']
    table['interactions'] = table['interactions'].astype(int)

    return table.reset_index()


def draw(config, scores, output_stem):
    """Two columns on a shared row grid: the host with the most parasites on the left,
    the others stacked on the right, so every row is the same height."""
    figure_style.apply()
    taxon_colors = config['parasite_groups']
    pairs = scores[['host', 'parasite', 'group', 'niche']].drop_duplicates()
    hosts = list(dict.fromkeys(pairs['host']))
    sizes = {h: max((pairs['host'] == h).sum(), MIN_COLUMN) for h in hosts}
    first = max(hosts, key=sizes.get)
    others = [h for h in hosts if h != first]
    gap = 2
    n_rows = max(sizes[first], sum(sizes[h] for h in others) + gap * (len(others) - 1))
    # the legends take this much under the axes; the rows shrink only when a page would
    # not hold them at the shared row height
    foot = 1.55
    row = min(figure_style.ROW, (figure_style.MAX_HEIGHT - foot - 0.35) / n_rows)
    fig = plt.figure(figsize=(figure_style.WIDTH, row * n_rows + foot + 0.35))
    grid = fig.add_gridspec(n_rows, 2, wspace=1.45,
                            top=1 - 0.35 / fig.get_figheight(), bottom=foot / fig.get_figheight(),
                            left=0.19, right=0.94)
    axes = {first: fig.add_subplot(grid[:sizes[first], 0])}
    start = 0
    for host in others:
        axes[host] = fig.add_subplot(grid[start:start + sizes[host], 1])
        start += sizes[host] + gap

    for host, ax in axes.items():
        rows = pairs[pairs['host'] == host].reset_index(drop=True)
        # a box per parasite, in the colour of its group
        for i, row in rows.iterrows():
            values = scores.loc[(scores['host'] == host) & (scores['parasite'] == row['parasite']), 'score']
            color = taxon_colors[row['group']]
            ax.boxplot(values, positions=[i], vert=False, widths=0.62, whis=1.5,
                       showfliers=False, patch_artist=True,
                       boxprops=dict(facecolor=color, edgecolor='#333333', linewidth=0.6),
                       medianprops=dict(color='#111111', linewidth=1),
                       whiskerprops=dict(color='#333333', linewidth=0.6),
                       capprops=dict(color='#333333', linewidth=0.6))
        # the group and niche strips beside the names
        for i, row in rows.iterrows():
            for offset, color in ((0, taxon_colors[row['group']]), (1, NICHE_COLORS[row['niche']])):
                ax.add_patch(Rectangle((-0.8 - 0.035 * offset, i - 0.4), 0.022, 0.8, color=color,
                                       linewidth=0, clip_on=False, transform=ax.get_yaxis_transform()))
        ax.set_title(host, fontweight='bold', loc='left', pad=3)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([figure_style.short_name(p) for p in rows['parasite']], style='italic')
        ax.set_ylim(sizes[host] - 0.4, -0.6)
        ax.set_xlim(SCORE_FLOOR - 0.02, 1)
        ax.tick_params(axis='y', length=0)
        ax.grid(axis='x', color='#e6e6e6', linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ('top', 'right', 'left'):
            ax.spines[side].set_visible(False)
        ax.spines['bottom'].set_color('#999999')
    for host in (first, others[-1]):
        axes[host].set_xlabel('Confidence score')

    figure_style.stack_legends(fig, left=0.19, blocks=[
        ('Taxonomic group', [Patch(color=c, label=g) for g, c in taxon_colors.items()
                             if g in set(pairs['group'])]),
        ('Niche', [Patch(color=c, label=niche) for niche, c in NICHE_COLORS.items()])])
    for extension in ('svg', 'png'):
        fig.savefig(f'{output_stem}.{extension}', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--table', default='paper/tables/confidence.csv')
    parser.add_argument('--figure', default='paper/figures/confidence', help='output path without extension')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    scores = collect(config, args.data_dir)
    table = summarise(scores)
    table.to_csv(args.table, index=False, float_format='%.3f')
    print(f"{len(scores):,} interactions over {len(table)} host-parasite pairs, median score "
          f"{scores['score'].median():.3f}; wrote {args.table}")
    draw(config, scores, args.figure)
    print(f"Wrote {args.figure}.svg and .png")
