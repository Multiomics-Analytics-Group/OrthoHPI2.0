'''
One host-parasite pair from parasite protein to host cell: the predicted interactions of
a parasite with a host as a matrix of parasite protein families against host proteins,
and under it, sharing the columns, the cell types of the infected tissues expressing
each host protein. A parasite family is the orthology group of the parasite proteins,
which the predictions are made at the level of, named by the most common description
of its members. Writes paper/figures/case_study.pdf/.svg and the two matrices to
paper/tables/case_study_links.csv and paper/tables/case_study_cell_types.csv.
'''
import argparse
import os
import re
import sys
from collections import Counter

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils
from scripts import figure_style

DOT_SCALE = ['#deebf7', '#9ecae1', '#6baed6', '#3182bd', '#08519c']
HEAT_SCALE = ['#ffffff', '#fee6ce', '#fdae6b', '#e6550d', '#a63603']
TISSUE_COLORS = ['#0072B2', '#009E73', '#E69F00', '#CC79A7', '#56B4E9', '#D55E00']
# what STRING's descriptions say when they say nothing
UNNAMED = {'uncharacterized protein', 'hypothetical protein', ''}
# most characters of a family name before it is cut to an ellipsis
NAME_CHARS = 36
# a cell type expressing a smaller share of the host proteins is an empty row
MIN_EXPRESSED = 0.1


def family_name(descriptions):
    '''The description most of the family's proteins carry, trimmed to what it names.'''
    counts = Counter(re.sub(r'[.\s]+$', '', str(d)).strip() for d in descriptions)
    named = [(n, d) for d, n in counts.most_common() if d.lower() not in UNNAMED]
    if not named:
        return 'uncharacterized'
    name = named[0][1]
    name = re.sub(r'\s*(domain[- ]containing protein|domain protein|family protein|protein)$',
                  '', name, flags=re.IGNORECASE)
    name = name[0].upper() + name[1:]

    return name if len(name) <= NAME_CHARS else name[:NAME_CHARS - 1].rstrip() + '\u2026'


def build(config, data_dir, parasite, host, min_targets, cells_per_tissue):
    '''The link matrix (families x host proteins) and the cell-type matrix (cell types x
    host proteins) of one pair, both keyed by host protein id. Families reaching fewer
    than min_targets host proteins are left out, and each tissue keeps the
    cells_per_tissue cell types expressing the most of the host proteins.'''
    predictions = pd.read_parquet(os.path.join(data_dir, 'predictions.parquet'))
    predictions['weight'] = predictions['weight'].astype(float)
    annotations = pd.read_parquet(os.path.join(data_dir, 'protein_annotations.parquet'))
    descriptions = annotations.set_index('protein')['description']
    edges = predictions[(predictions['taxid1'] == str(parasite))
                        & (predictions['taxid2'] == str(host))].copy()
    if edges.empty:
        raise SystemExit(f'no predictions for parasite {parasite} in host {host}')

    edges['description'] = edges['source'].map(descriptions)
    families = edges.groupby('group1').agg(
        proteins=('source', 'nunique'), targets=('target', 'nunique'),
        name=('description', family_name)).sort_values(['targets', 'proteins'], ascending=False)
    families = families[families['targets'] >= min_targets]
    # an unnamed family keeps its group id, so two of them are told apart
    families['label'] = [(f'{row.name} ({row.proteins})' if row.proteins > 1 else row.name)
                         + (f' {row.Index}' if row.name == 'uncharacterized' else '')
                         for row in families.itertuples()]
    edges = edges[edges['group1'].isin(families.index)]
    links = edges.groupby(['group1', 'target'])['weight'].max().unstack()
    # columns by host family, so paralogues sit together, families by their reach
    symbols = edges.drop_duplicates('target').set_index('target')['target_name'].astype(str).str.upper()
    host_groups = edges.drop_duplicates('target').set_index('target')['group2']
    column_order = sorted(links.columns, key=lambda t: (-int(links[t].notna().sum()), host_groups[t], symbols[t]))
    links = links.reindex(index=families.index, columns=column_order)

    tissue_labels = [config['tissues'][t] for t in config['parasites'][int(parasite)]['tissues']]
    cells = pd.read_parquet(os.path.join(data_dir, 'tissues_cell_types.parquet'))
    cells = cells[cells['Gene'].isin(links.columns) & cells['Tissue'].isin(tissue_labels)]
    cells = cells.dropna(subset=['Cell type'])
    expression = (cells.groupby(['Tissue', 'Cell type', 'Gene'])['nTPM'].max().unstack()
                  .reindex(columns=links.columns))
    # tissues in the parasite's order; cell types by how many of the targets they express
    expressed = (expression > 0).sum(axis=1)
    expression = expression[expressed >= MIN_EXPRESSED * len(links.columns)]
    expression = expression.loc[sorted(expression.index,
                                       key=lambda i: (tissue_labels.index(i[0]), -expressed[i], i[1]))]
    expression = expression.groupby(level=0, sort=False).head(cells_per_tissue)

    return links, families, symbols, expression, tissue_labels


def draw(links, families, symbols, expression, tissue_labels, title, output_stem):
    figure_style.apply()
    n_fam, n_cols, n_cells = len(links), links.shape[1], len(expression)
    fig_height = 0.15 * (n_fam + n_cells) + 2.3
    fig, (top, bottom) = plt.subplots(
        2, 1, figsize=(figure_style.WIDTH, fig_height), sharex=True,
        gridspec_kw={'height_ratios': [n_fam, n_cells], 'hspace': 0.04})

    # A: a dot per predicted link, shaded by confidence
    rows, cols = np.where(links.notna().to_numpy())
    scores = links.to_numpy()[rows, cols]
    cmap = LinearSegmentedColormap.from_list('score', DOT_SCALE)
    dots = top.scatter(cols, rows, c=scores, cmap=cmap, vmin=0.35, vmax=1.0, s=22,
                       linewidths=0, zorder=3)
    top.set_yticks(range(n_fam))
    top.set_yticklabels(families['label'])
    top.set_ylim(n_fam - 0.5, -0.5)
    top.set_title('A  Parasite families and the host proteins they are predicted to bind',
                  loc='left', pad=3)

    # B: the cell types of the infected tissues expressing each host protein
    values = expression.to_numpy(dtype=float)
    values = np.where(values > 0, values, np.nan)
    heat_cmap = LinearSegmentedColormap.from_list('ntpm', HEAT_SCALE)
    heat_cmap.set_bad('#ffffff')
    floor = max(np.nanmin(values), 0.1)
    heat = bottom.imshow(np.ma.masked_invalid(values), cmap=heat_cmap, aspect='auto',
                         norm=LogNorm(vmin=floor, vmax=np.nanmax(values)), interpolation='nearest')
    bottom.set_yticks(range(n_cells))
    bottom.set_yticklabels([cell for _, cell in expression.index])
    bottom.set_ylim(n_cells - 0.5, -0.5)
    bottom.set_xticks(range(n_cols))
    bottom.set_xticklabels([symbols[t] for t in links.columns], rotation=90)
    bottom.set_title('B  Cell types of the infected tissues expressing them', loc='left', pad=3)
    # a tissue strip down the side of B, one colour per tissue of the parasite
    tissue_color = dict(zip(tissue_labels, TISSUE_COLORS))
    for i, (tissue, _) in enumerate(expression.index):
        bottom.add_patch(Rectangle((-1.5, i - 0.5), 0.7, 1, color=tissue_color[tissue],
                                   linewidth=0, clip_on=False))
    bottom.set_xlim(-0.5, n_cols - 0.5)

    for ax in (top, bottom):
        ax.tick_params(length=0, pad=2)
        for side in ax.spines.values():
            side.set_visible(False)
        ax.set_xticks(np.arange(-0.5, n_cols), minor=True)
        ax.tick_params(which='minor', length=0)
    # the cell-type names sit left of the tissue strip
    bottom.tick_params(axis='y', pad=11)
    top.grid(color='#eeeeee', linewidth=0.5)
    top.set_axisbelow(True)
    bottom.grid(which='minor', axis='x', color='#ffffff', linewidth=0.4)

    # scales beside the panels, the tissue legend under them
    score_bar = fig.colorbar(dots, ax=top, fraction=0.02, pad=0.01)
    score_bar.set_label('confidence score')
    score_bar.ax.tick_params(length=2)
    score_bar.outline.set_visible(False)
    heat_bar = fig.colorbar(heat, ax=bottom, fraction=0.02, pad=0.01)
    heat_bar.set_label('nTPM')
    heat_bar.ax.tick_params(length=2)
    heat_bar.outline.set_visible(False)
    handles = [Patch(color=tissue_color[t], label=t) for t in tissue_labels if t in {i[0] for i in expression.index}]
    figure_style.stack_legends(fig, left=0.36, blocks=[('Tissue', handles)], ncol=len(handles))
    fig.suptitle(title, x=0.02, ha='left', y=0.995, style='italic')
    fig.subplots_adjust(left=0.36, right=0.92, top=1 - 0.45 / fig_height, bottom=1.1 / fig_height)
    for extension in ('pdf', 'svg'):
        fig.savefig(f'{output_stem}.{extension}')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--parasite', type=int, default=51031,
                        help='parasite taxid (default: Necator americanus)')
    parser.add_argument('--host', type=int, default=9606, help='host taxid (default: human)')
    parser.add_argument('--min-targets', type=int, default=3,
                        help='families reaching fewer host proteins are left out')
    parser.add_argument('--cells-per-tissue', type=int, default=8,
                        help='the cell types of each tissue kept: those expressing the most host proteins')
    parser.add_argument('--links', default='paper/tables/case_study_links.csv')
    parser.add_argument('--cell-types', default='paper/tables/case_study_cell_types.csv')
    parser.add_argument('--figure', default='paper/figures/case_study',
                        help='output path without extension')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    links, families, symbols, expression, tissues = build(
        config, args.data_dir, args.parasite, args.host, args.min_targets, args.cells_per_tissue)
    (links.set_axis(families['label'], axis=0).set_axis([symbols[t] for t in links.columns], axis=1)
     .to_csv(args.links))
    expression.set_axis([symbols[t] for t in expression.columns], axis=1).to_csv(args.cell_types)
    print(f"{len(families)} parasite families, {links.shape[1]} host proteins, "
          f"{len(expression)} cell types; wrote {args.links} and {args.cell_types}")
    title = (f"{config['parasites'][args.parasite]['label']} in "
             f"{config['hosts'][args.host]['label'].split(' (')[0]}")
    draw(links, families, symbols, expression, tissues, title, args.figure)
    print(f"Wrote {args.figure}.pdf and .svg")
