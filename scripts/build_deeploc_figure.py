'''
The DeepLoc 2 probabilities the localization filters kept proteins on, as on the app's
home page, in two figures. Host proteins: one panel per class, a box per host over the
proteins called for that class; the cytoplasm and nucleus panels hold only the host
proteins of the intracellular parasites, the only ones the filter reads those classes
for. Parasite proteins: a box per parasite over the extracellular probability of its
proteins called extracellular, and over the cell-membrane probability of those called
cell membrane for the unicellular parasites, the secretome filter keeping only the
secreted proteins of a multicellular one. Every prediction is counted and each protein
once. Writes paper/tables/deeploc.csv (one row per box) and
paper/figures/deeploc_host.svg/.png and deeploc_parasite.svg/.png.
'''
import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch, Rectangle
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils
from pipeline import main
from scripts import figure_style

NICHE_COLORS = {'extracellular': '#c7c7c7', 'intracellular': '#3d3d3d'}
# how DeepLoc's class names are spelt in the localisations table
COLUMN_OF = {'Extracellular': 'extracellular', 'Cell membrane': 'cell_membrane',
             'Cytoplasm': 'cytoplasm', 'Nucleus': 'nucleus'}
SURFACE = ['Extracellular', 'Cell membrane']
# the outline colours of the classes in the app
CLASS_COLORS = {'Extracellular': '#3690c0', 'Cell membrane': '#045a8d',
                'Cytoplasm': '#e6550d', 'Nucleus': '#a63603'}
# the probability axes start a little under the lowest cut-off
FLOOR = min(main.DEEPLOC_CUTOFFS.values()) - 0.03


def load(config, data_dir):
    '''One row per predicted protein, host and parasite, with the side it is on, the
    host, the parasite's group and niche, and its DeepLoc probabilities.'''
    localisations = pd.read_parquet(os.path.join(data_dir, 'deeploc_localisations.parquet'))
    predictions = pd.read_parquet(os.path.join(data_dir, 'predictions.parquet'))
    hosts = {str(t): h['label'].split('(')[1].rstrip(')') for t, h in config['hosts'].items()}
    parasites = {str(t): p for t, p in config['parasites'].items()}

    frames = []
    for side in ('source', 'target'):
        frame = predictions[['taxid1', 'taxid2', side]].drop_duplicates().rename(
            columns={side: 'protein'})
        frames.append(frame.assign(side='parasite' if side == 'source' else 'host'))
    proteins = pd.concat(frames, ignore_index=True)
    proteins['host'] = proteins['taxid2'].map(hosts)
    proteins['parasite'] = proteins['taxid1'].map(lambda t: parasites[t]['label'])
    proteins['group'] = proteins['taxid1'].map(lambda t: parasites[t]['group'])
    proteins['niche'] = proteins['taxid1'].map(lambda t: parasites[t]['niche'])
    proteins['multicellular'] = proteins['taxid1'].map(lambda t: parasites[t].get('multicellular', False))

    return proteins.merge(localisations, on='protein', how='inner')


def host_scores(proteins, config):
    '''The probability of each host protein for every class it was called for and some
    parasite reaches it in, once per host and class; hosts in config order.'''
    hosts = proteins[proteins['side'] == 'host']
    rows = []
    for surface_class, cutoff in main.DEEPLOC_CUTOFFS.items():
        column = COLUMN_OF[surface_class]
        allowed = hosts if surface_class in SURFACE else hosts[hosts['niche'] == 'intracellular']
        called = allowed[allowed[column] > cutoff].drop_duplicates(['host', 'protein'])
        rows.append(called.assign(localization=surface_class, probability=called[column]))
    order = [h['label'].split('(')[1].rstrip(')') for h in config['hosts'].values()]
    scores = pd.concat(rows, ignore_index=True)
    scores['host'] = pd.Categorical(scores['host'], categories=order, ordered=True)

    return scores[['localization', 'host', 'protein', 'probability']].sort_values(
        'host', kind='stable')


def parasite_scores(proteins, config):
    '''The extracellular probability of each parasite protein called extracellular, and
    the cell-membrane probability of each unicellular parasite protein called cell
    membrane, once per parasite; parasites by group, niche then name.'''
    parasites = proteins[proteins['side'] == 'parasite'].drop_duplicates(['parasite', 'protein'])
    rows = []
    for surface_class in SURFACE:
        column = COLUMN_OF[surface_class]
        allowed = parasites if surface_class == 'Extracellular' else parasites[~parasites['multicellular']]
        called = allowed[allowed[column] > main.DEEPLOC_CUTOFFS[surface_class]]
        rows.append(called.assign(localization=surface_class, probability=called[column]))
    scores = pd.concat(rows, ignore_index=True)
    groups = list(config['parasite_groups'])
    scores['rank'] = list(zip(scores['group'].map(groups.index),
                              scores['niche'].map(list(NICHE_COLORS).index), scores['parasite']))

    return scores.sort_values('rank')[['localization', 'parasite', 'group', 'niche',
                                       'protein', 'probability']]


def summarise(host, parasite):
    '''One row per box: the proteins in it and the quartiles of their probabilities.'''
    tables = []
    for side, scores, by in (('host', host, 'host'), ('parasite', parasite, 'parasite')):
        boxes = scores.groupby(['localization', by], sort=False, observed=True)['probability']
        table = boxes.describe()[['count', 'min', '25%', '50%', '75%', 'max']]
        table.columns = ['proteins', 'min', 'q1', 'median', 'q3', 'max']
        table['proteins'] = table['proteins'].astype(int)
        tables.append(table.reset_index().rename(columns={by: 'species'}).assign(side=side))

    return pd.concat(tables, ignore_index=True)[['side', 'localization', 'species', 'proteins',
                                                 'min', 'q1', 'median', 'q3', 'max']]


def draw_box(ax, values, position, color, points, rng):
    '''A box over `values` with the proteins as jittered points behind it.'''
    ax.scatter(values, position + rng.uniform(-0.28, 0.28, len(values)), s=points,
               color=color, alpha=0.35, linewidths=0, zorder=1)
    ax.boxplot(values, positions=[position], vert=False, widths=0.62, whis=1.5,
               showfliers=False, patch_artist=True, zorder=2,
               boxprops=dict(facecolor='none', edgecolor='#222222', linewidth=0.7),
               medianprops=dict(color='#111111', linewidth=1.1),
               whiskerprops=dict(color='#222222', linewidth=0.7),
               capprops=dict(color='#222222', linewidth=0.7))


def style_axis(ax, cutoff, color, n_rows):
    ax.axvline(cutoff, color=color, linestyle=':', linewidth=1, zorder=0)
    ax.set_xlim(FLOOR, 1.005)
    ax.set_ylim(n_rows - 0.4, -0.6)
    ax.tick_params(axis='y', length=0)
    ax.grid(axis='x', color='#e6e6e6', linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_color('#999999')


def draw_hosts(scores, output_stem):
    '''Four panels, one per class, a row per host.'''
    figure_style.apply()
    hosts = list(scores['host'].cat.categories)
    fig, axes = plt.subplots(2, 2, figsize=(figure_style.WIDTH, 3.4), sharey=True)
    rng = np.random.default_rng(0)
    for letter, ax, (surface_class, cutoff) in zip('ABCD', axes.flat, main.DEEPLOC_CUTOFFS.items()):
        color = CLASS_COLORS[surface_class]
        for i, host in enumerate(hosts):
            values = scores.loc[(scores['localization'] == surface_class) & (scores['host'] == host),
                                'probability']
            if len(values):
                draw_box(ax, values, i, color, 2, rng)
            else:
                ax.text(FLOOR + 0.01, i, 'no proteins', va='center', color='#888888')
        style_axis(ax, cutoff, color, len(hosts))
        ax.set_yticks(range(len(hosts)))
        ax.set_yticklabels(hosts)
        ax.set_title(f'({letter}) {surface_class}', fontweight='bold', loc='left', pad=3)
        ax.set_xlabel(f'P({surface_class.lower()})')
    fig.tight_layout(h_pad=1.2, w_pad=2)
    for extension in ('svg', 'png'):
        fig.savefig(f'{output_stem}.{extension}', dpi=300)
    plt.close(fig)


def draw_parasites(config, scores, output_stem):
    '''Two panels side by side: every parasite's extracellular proteins on the left, the
    unicellular parasites' membrane proteins on the right, on a shared row grid.'''
    figure_style.apply()
    taxon_colors = config['parasite_groups']
    panels = {c: scores[scores['localization'] == c] for c in SURFACE}
    rows = {c: list(dict.fromkeys(panels[c]['parasite'])) for c in SURFACE}
    n_rows = max(len(r) for r in rows.values())
    foot = 1.55
    row = min(figure_style.ROW, (figure_style.MAX_HEIGHT - foot - 0.35) / n_rows)
    fig = plt.figure(figsize=(figure_style.WIDTH, row * n_rows + foot + 0.35))
    grid = fig.add_gridspec(n_rows, 2, wspace=1.45,
                            top=1 - 0.35 / fig.get_figheight(), bottom=foot / fig.get_figheight(),
                            left=0.19, right=0.94)
    rng = np.random.default_rng(0)
    for column, (letter, surface_class) in enumerate(zip('AB', SURFACE)):
        names = rows[surface_class]
        ax = fig.add_subplot(grid[:len(names), column])
        panel = panels[surface_class]
        for i, parasite in enumerate(names):
            box = panel[panel['parasite'] == parasite]
            first = box.iloc[0]
            draw_box(ax, box['probability'], i, taxon_colors[first['group']], 2, rng)
            # the group and niche strips beside the names
            for offset, color in ((0, taxon_colors[first['group']]), (1, NICHE_COLORS[first['niche']])):
                ax.add_patch(Rectangle((-0.8 - 0.035 * offset, i - 0.4), 0.022, 0.8, color=color,
                                       linewidth=0, clip_on=False, transform=ax.get_yaxis_transform()))
        style_axis(ax, main.DEEPLOC_CUTOFFS[surface_class], '#777777', len(names))
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels([figure_style.short_name(p) for p in names], style='italic')
        ax.set_title(f'({letter}) {surface_class}', fontweight='bold', loc='left', pad=3)
        ax.set_xlabel(f'P({surface_class.lower()})')

    figure_style.stack_legends(fig, left=0.19, blocks=[
        ('Taxonomic group', [Patch(color=c, label=g) for g, c in taxon_colors.items()
                             if g in set(scores['group'])]),
        ('Niche', [Patch(color=c, label=niche) for niche, c in NICHE_COLORS.items()])])
    for extension in ('svg', 'png'):
        fig.savefig(f'{output_stem}.{extension}', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--table', default='paper/tables/deeploc.csv')
    parser.add_argument('--figures', default='paper/figures/deeploc',
                        help='output path without the _host/_parasite suffix and extension')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    proteins = load(config, args.data_dir)
    host = host_scores(proteins, config)
    parasite = parasite_scores(proteins, config)
    table = summarise(host, parasite)
    table.to_csv(args.table, index=False, float_format='%.3f')
    print(f"{len(host):,} host and {len(parasite):,} parasite protein-class pairs; wrote {args.table}")
    draw_hosts(host, f'{args.figures}_host')
    draw_parasites(config, parasite, f'{args.figures}_parasite')
    print(f"Wrote {args.figures}_host and {args.figures}_parasite .svg and .png")
