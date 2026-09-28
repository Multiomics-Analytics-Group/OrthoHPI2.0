'''
The predicted networks at a glance: one block per host, as tall as the number of
parasites infecting it, a bar per parasite, the blocks sharing the axis so a bar is
comparable across hosts. Each bar is split by where DeepLoc puts the host protein -- on
the surface (extracellular or cell membrane), inside the cell (cytoplasm or nucleus) or
called for both -- which is what the parasite's niche opened to it; the taxonomic group
and niche are the strips beside the names. Writes paper/tables/networks.csv (one row per
host-parasite pair) and paper/figures/networks.pdf/.svg/.png.
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
from pipeline import main
from scripts import figure_style

NICHE_COLORS = {'extracellular': '#c7c7c7', 'intracellular': '#3d3d3d'}
PLACES = ['surface', 'both', 'interior']
PLACE_LABELS = {'surface': 'surface (extracellular, cell membrane)',
                'both': 'surface and interior', 'interior': 'interior (cytoplasm, nucleus)'}
PLACE_COLORS = {'surface': '#9ecae1', 'both': '#4292c6', 'interior': '#08519c'}
# how DeepLoc's class names are spelt in the localisations table
COLUMN_OF = {'Extracellular': 'extracellular', 'Cell membrane': 'cell_membrane',
             'Cytoplasm': 'cytoplasm', 'Nucleus': 'nucleus'}
# a host with fewer parasites still gets a block this many bars tall
MIN_COLUMN = 3


def host_places(data_dir):
    '''{host protein: surface / interior / both}, on the pipeline's own cut-offs.'''
    localisations = pd.read_parquet(os.path.join(data_dir, 'deeploc_localisations.parquet'))
    called = {COLUMN_OF[c]: localisations[COLUMN_OF[c]] > cutoff
              for c, cutoff in main.DEEPLOC_CUTOFFS.items()}
    surface = called['extracellular'] | called['cell_membrane']
    interior = called['cytoplasm'] | called['nucleus']
    place = pd.Series('surface', index=localisations['protein'].values)
    place[(surface & interior).values] = 'both'
    place[(~surface & interior).values] = 'interior'

    return place


def count(config, data_dir):
    '''One row per host-parasite pair: interactions, by host-protein place, and the
    proteins on either side; hosts in config order, parasites by group then name.'''
    parasites = config['parasites']
    hosts = {str(t): h['label'].split('(')[1].rstrip(')') for t, h in config['hosts'].items()}
    predictions = pd.read_parquet(os.path.join(data_dir, 'predictions.parquet'))
    predictions['place'] = predictions['target'].map(host_places(data_dir))
    order = list(config['parasite_groups'])

    rows = []
    for host, host_name in hosts.items():
        for taxid, parasite in sorted(parasites.items(),
                                      key=lambda item: (order.index(item[1]['group']), item[1]['label'])):
            if int(host) not in parasite['hosts']:
                continue
            edges = predictions[(predictions['taxid1'] == str(taxid)) & (predictions['taxid2'] == host)]
            by_place = edges.groupby('place').size()
            rows.append({'host': host_name, 'parasite': parasite['label'], 'group': parasite['group'],
                         'niche': parasite['niche'], 'interactions': len(edges),
                         'parasite proteins': edges['source'].nunique(),
                         'host proteins': edges['target'].nunique(),
                         **{place: int(by_place.get(place, 0)) for place in PLACES}})

    return pd.DataFrame(rows)


def draw(config, table, output_stem):
    """Two columns on a shared row grid: the host with the most parasites on the left,
    the others stacked on the right, so every row is the same height."""
    figure_style.apply()
    taxon_colors = config['parasite_groups']
    hosts = list(dict.fromkeys(table['host']))
    sizes = {h: max((table['host'] == h).sum(), MIN_COLUMN) for h in hosts}
    first = max(hosts, key=sizes.get)
    others = [h for h in hosts if h != first]
    gap = 2
    n_rows = max(sizes[first], sum(sizes[h] for h in others) + gap * (len(others) - 1))
    # the legends take this much under the axes; the rows shrink only when a page would
    # not hold them at the shared row height
    foot = 1.95
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
        rows = table[table['host'] == host].reset_index(drop=True)
        y = range(len(rows))
        # a bar per parasite, split by where the host protein is
        left = pd.Series(0.0, index=rows.index)
        for place in PLACES:
            ax.barh(y, rows[place], left=left, color=PLACE_COLORS[place], height=0.72)
            left += rows[place]
        for i, row in rows.iterrows():
            ax.text(row['interactions'] + 20, i, f"{row['interactions']:,}",
                    va='center', color='#555555')
        # the group and niche strips beside the names
        for i, row in rows.iterrows():
            for offset, color in ((0, taxon_colors[row['group']]), (1, NICHE_COLORS[row['niche']])):
                ax.add_patch(Rectangle((-0.8 - 0.035 * offset, i - 0.4), 0.022, 0.8, color=color,
                                       linewidth=0, clip_on=False, transform=ax.get_yaxis_transform()))
        ax.set_title(host, fontweight='bold', loc='left', pad=3)
        ax.set_yticks(list(y))
        ax.set_yticklabels([figure_style.short_name(p) for p in rows['parasite']], style='italic')
        ax.set_ylim(sizes[host] - 0.4, -0.6)
        ax.set_xlim(0, table['interactions'].max() * 1.2)
        ax.tick_params(length=0)
        ax.grid(axis='x', color='#e6e6e6', linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ('top', 'right', 'left'):
            ax.spines[side].set_visible(False)
        ax.spines['bottom'].set_color('#999999')
    for host in (first, others[-1]):
        axes[host].set_xlabel('Predicted interactions')

    figure_style.stack_legends(fig, left=0.19, blocks=[
        ('Host-protein location', [Patch(color=PLACE_COLORS[p], label=PLACE_LABELS[p]) for p in PLACES], 2),
        ('Taxonomic group', [Patch(color=c, label=g) for g, c in taxon_colors.items()
                             if g in set(table['group'])]),
        ('Niche', [Patch(color=c, label=niche) for niche, c in NICHE_COLORS.items()])])
    for extension in ('pdf', 'svg', 'png'):
        fig.savefig(f'{output_stem}.{extension}', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--table', default='paper/tables/networks.csv')
    parser.add_argument('--figure', default='paper/figures/networks', help='output path without extension')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    table = count(config, args.data_dir)
    table.to_csv(args.table, index=False)
    print(f"{table['interactions'].sum():,} interactions over {len(table)} host-parasite pairs; "
          f"wrote {args.table}")
    draw(config, table, args.figure)
    print(f"Wrote {args.figure}.pdf, .svg and .png")
