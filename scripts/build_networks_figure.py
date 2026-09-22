'''
The predicted networks at a glance: one block per host, as tall as the number of
parasites infecting it, a bar per parasite, the hosts sharing the axis so a bar is
comparable across them. Each bar is split by where DeepLoc puts the host protein -- on
the surface (extracellular or cell membrane), inside the cell (cytoplasm or nucleus) or
called for both -- which is what the parasite's niche opened to it; the taxonomic group
and niche are the strips beside the names. Writes paper/tables/networks.csv (one row per
host-parasite pair) and paper/figures/networks.pdf/.svg.
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
    matplotlib.rcParams['font.family'] = 'sans-serif'
    taxon_colors = config['parasite_groups']
    hosts = list(dict.fromkeys(table['host']))
    heights = [max((table['host'] == h).sum(), MIN_COLUMN) for h in hosts]
    fig, axes = plt.subplots(len(hosts), 1, figsize=(7.5, 0.15 * sum(heights) + 1.6), sharex=True,
                             gridspec_kw={'height_ratios': heights, 'hspace': 0.12})

    for ax, host in zip(axes, hosts):
        rows = table[table['host'] == host].reset_index(drop=True)
        y = range(len(rows))
        # a bar per parasite, split by where the host protein is
        left = pd.Series(0.0, index=rows.index)
        for place in PLACES:
            ax.barh(y, rows[place], left=left, color=PLACE_COLORS[place], height=0.72)
            left += rows[place]
        for i, row in rows.iterrows():
            ax.text(row['interactions'] + 12, i, f"{row['interactions']:,}", fontsize=5.5,
                    va='center', color='#555555')
        # the group and niche strips beside the names
        for i, row in rows.iterrows():
            for offset, color in ((0, taxon_colors[row['group']]), (1, NICHE_COLORS[row['niche']])):
                ax.add_patch(Rectangle((-0.29 - 0.012 * offset, i - 0.4), 0.008, 0.8, color=color,
                                       linewidth=0, clip_on=False, transform=ax.get_yaxis_transform()))
        ax.set_title(host, fontsize=8, fontweight='bold', loc='left', pad=3)
        ax.set_yticks(list(y))
        ax.set_yticklabels(rows['parasite'], fontsize=6, style='italic')
        ax.set_ylim(max(len(rows), MIN_COLUMN) - 0.4, -0.6)
        ax.tick_params(length=0)
        ax.tick_params(axis='x', labelsize=6, labelbottom=host == hosts[-1])
        ax.grid(axis='x', color='#e6e6e6', linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ('top', 'right', 'left'):
            ax.spines[side].set_visible(False)
        ax.spines['bottom'].set_color('#999999')
    axes[0].set_xlim(0, table['interactions'].max() * 1.08)
    axes[-1].set_xlabel('Predicted interactions', fontsize=7)

    handles = [Patch(color=PLACE_COLORS[p], label=PLACE_LABELS[p]) for p in PLACES]
    handles += [Patch(color=c, label=g) for g, c in taxon_colors.items() if g in set(table['group'])]
    handles += [Patch(color=c, label=niche) for niche, c in NICHE_COLORS.items()]
    fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=6, frameon=False,
               handlelength=0.9, bbox_to_anchor=(0.5, 0.0), columnspacing=1.2)
    fig.subplots_adjust(left=0.27, right=0.97, top=1 - 0.3 / fig.get_figheight(),
                        bottom=0.85 / fig.get_figheight())
    for extension in ('pdf', 'svg'):
        fig.savefig(f'{output_stem}.{extension}')
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
    print(f"Wrote {args.figure}.pdf and .svg")
