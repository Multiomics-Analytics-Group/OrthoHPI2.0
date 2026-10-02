'''
How broadly the predicted human interactors are expressed across the cell types of the
tissues they are predicted in, against every HPA gene detected in the same tissues.

For each (human target, infected tissue) pair, the breadth is the share of the tissue's HPA
cell types where the target's nTPM is above the cutoff; pairs the target is not detected in
at all, or whose tissue HPA does not cover, are left out, and the background is held to the
same rule.

Usage: .venv/bin/python scripts/plot_cell_type_breadth.py [--ntpm 1] [--output FILE]
'''
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
import figure_style
import utils
from plot_hpa_ntpm_distribution import read_hpa_with_zeros

PREDICTED_COLOR = '#2a78d6'
BACKGROUND_COLOR = '#a8a79f'
INK_SECONDARY = '#52514e'
GRID = '#dcdbd5'
BINS = [0, 0.25, 0.5, 0.75, 0.9999, 1.0]
BIN_LABELS = ['≤25%', '≤50%', '≤75%', '<100%', 'All']


def predicted_pairs(config):
    '''The (human target, tissue) pairs of the predictions, over each parasite's tissues.'''
    tissue_names = {bto: name.lower() for bto, name in config['tissues'].items()}
    predictions = pd.read_parquet('data/predictions.parquet')
    predictions = predictions[predictions['taxid2'] == '9606']
    rows = [(target, tissue_names[bto])
            for taxid, group in predictions.groupby('taxid1')
            for bto in config['parasites'][int(taxid)]['tissues']
            for target in group['target'].unique()]

    return pd.DataFrame(rows, columns=['Gene', 'Tissue']).drop_duplicates()


def breadth(hpa, cutoff):
    '''Per (gene, tissue): the share of the tissue's cell types expressing the gene.'''
    cell_types = hpa.groupby('Tissue')['Cell type'].nunique()
    data = (hpa.assign(on=hpa['nTPM'] > cutoff)
            .groupby(['Gene', 'Tissue'])['on'].sum().rename('n_on').reset_index())
    data['n_cell'] = data['Tissue'].map(cell_types)
    data['share'] = data['n_on'] / data['n_cell']

    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--ntpm', type=float, default=1.0)
    parser.add_argument('--output', default='paper/figures/cell_type_breadth.pdf')
    args = parser.parse_args()

    config = utils.read_config(filepath='config.yml')
    per_gene = breadth(read_hpa_with_zeros('config.yml'), args.ntpm)
    per_gene = per_gene[per_gene['n_on'] > 0]
    pairs = predicted_pairs(config)
    predicted = pairs.merge(per_gene, on=['Gene', 'Tissue'])
    background = per_gene[per_gene['Tissue'].isin(predicted['Tissue'].unique())]

    for name, data in [('predicted', predicted), ('background', background)]:
        print(f'{name}: {len(data)} pairs, in all cell types {(data["share"] == 1).mean():.1%}, '
              f'in ≥75% {(data["share"] >= 0.75).mean():.1%}')

    figure_style.apply()
    fig, (left, right) = plt.subplots(1, 2, figsize=(figure_style.WIDTH, 3.4),
                                      gridspec_kw={'width_ratios': [1, 1.25], 'wspace': 0.55})

    # left: the breadth distribution, predicted beside background
    x = np.arange(len(BIN_LABELS))
    width = 0.38
    for offset, data, color, label in [(-width / 2, predicted, PREDICTED_COLOR, 'Predicted interactors'),
                                       (width / 2, background, BACKGROUND_COLOR, 'All detected HPA genes')]:
        shares = pd.cut(data['share'], BINS, labels=BIN_LABELS, include_lowest=True).value_counts(normalize=True)
        shares = shares.reindex(BIN_LABELS)
        left.bar(x + offset, shares * 100, width - 0.03, color=color, label=label, zorder=2)
        for xi, value in zip(x + offset, shares):
            left.text(xi, value * 100 + 1, f'{value:.0%}', ha='center', va='bottom', fontsize=6,
                      color=INK_SECONDARY)
    left.set_xticks(x, BIN_LABELS)
    left.set_xlabel(f"Share of the tissue's cell types\nexpressing the protein (nTPM > {args.ntpm:g})")
    left.set_ylabel('% of (protein, tissue) pairs')
    fig.legend(*left.get_legend_handles_labels(), frameon=False, loc='upper center', ncol=2,
               bbox_to_anchor=(0.5, 1.04))
    left.grid(axis='y', color=GRID, linewidth=0.6, zorder=0)

    # right: the share in all cell types, per tissue, ordered by how many cell types it has
    tissues = (predicted.groupby('Tissue').agg(n_cell=('n_cell', 'first'),
                                               predicted=('share', lambda s: (s == 1).mean()))
               .join(background.groupby('Tissue')['share'].apply(lambda s: (s == 1).mean()).rename('background'))
               .sort_values(['n_cell', 'predicted']))
    y = np.arange(len(tissues))
    right.hlines(y, tissues['background'] * 100, tissues['predicted'] * 100, color=GRID, linewidth=2, zorder=1)
    right.scatter(tissues['background'] * 100, y, s=30, color=BACKGROUND_COLOR, zorder=2,
                  edgecolor='white', linewidth=1)
    right.scatter(tissues['predicted'] * 100, y, s=30, color=PREDICTED_COLOR, zorder=3,
                  edgecolor='white', linewidth=1)
    right.set_yticks(y, [f'{tissue} ({n})' for tissue, n in tissues['n_cell'].items()])
    right.set_xlim(0, 100)
    right.set_xlabel('% of pairs expressed in every cell type')
    right.grid(axis='x', color=GRID, linewidth=0.6, zorder=0)

    for ax in (left, right):
        ax.spines[['top', 'right']].set_visible(False)
    left.set_title('A', loc='left', fontweight='bold')
    right.set_title('B', loc='left', fontweight='bold')
    right.set_title('(n): HPA cell types in the tissue', loc='right', color=INK_SECONDARY, fontsize=figure_style.FONT)

    fig.savefig(args.output, bbox_inches='tight', dpi=200)
    print(f'wrote {args.output}')


if __name__ == '__main__':
    main()
