'''
How alike the predicted interactomes of the parasites of one host are: the Jaccard
similarity of the host proteins each pair of parasites reaches, drawn as the heatmap of
the "Parasites of a host" page, with the taxonomic group and niche strips along both
axes and the parasites ordered by niche, then group, then name (or group first, with
--order-by group), so the two niches read as blocks. Writes
paper/figures/parasite_similarity.pdf/.svg, the matrix to
paper/tables/parasite_similarity.csv and, for the text, the median similarity within
and between genera, taxonomic groups, niches and organisations to
paper/tables/parasite_similarity_medians.csv.
'''
import argparse
import itertools
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils

NICHE_COLORS = {'extracellular': '#c7c7c7', 'intracellular': '#3d3d3d'}
NICHE_ORDER = ['extracellular', 'intracellular']
# the scale of the app's heatmap
SCALE = ['#ffffff', '#deebf7', '#9ecae1', '#6baed6', '#3182bd', '#08519c']
DIAGONAL_COLOR = '#d9d9d9'


def interactors(config, predictions, host, order_by):
    '''{parasite taxid: host proteins it reaches}, in the order the axes take: the chosen
    annotation first, so its values are contiguous, then the other, then the name.'''
    parasites = config['parasites']
    group_order = list(config['parasite_groups'])
    edges = predictions[predictions['taxid2'] == str(host)]
    targets = {int(t): set(rows['target']) for t, rows in edges.groupby('taxid1')}

    def rank(taxid):
        parasite = parasites[taxid]
        group, niche = group_order.index(parasite['group']), NICHE_ORDER.index(parasite['niche'])
        return ((niche, group, parasite['label']) if order_by == 'niche'
                else (group, niche, parasite['label']))

    return {taxid: targets[taxid] for taxid in sorted(targets, key=rank)}


def similarity_matrix(targets):
    taxids = list(targets)
    matrix = np.full((len(taxids), len(taxids)), np.nan)
    for i, a in enumerate(taxids):
        for j, b in enumerate(taxids):
            if i != j:
                matrix[i, j] = len(targets[a] & targets[b]) / len(targets[a] | targets[b])

    return pd.DataFrame(matrix, index=taxids, columns=taxids)


def medians(config, similarity):
    '''Median similarity of the pairs sharing, and not sharing, each annotation.'''
    parasites = config['parasites']
    rows = []
    for a, b in itertools.combinations(similarity.index, 2):
        pa, pb = parasites[a], parasites[b]
        rows.append({
            'similarity': similarity.loc[a, b],
            'genus': pa['label'].split()[0] == pb['label'].split()[0],
            'taxonomic group': pa['group'] == pb['group'],
            'niche': pa['niche'] == pb['niche'],
            'organisation': bool(pa.get('multicellular')) == bool(pb.get('multicellular'))})
    pairs = pd.DataFrame(rows)
    table = pd.DataFrame({
        'annotation': ['all'] + [c for c in pairs.columns if c != 'similarity'],
        'pairs': [len(pairs)] + [int(pairs[c].sum()) for c in pairs.columns if c != 'similarity'],
        'same': [pairs['similarity'].median()]
                + [pairs.loc[pairs[c], 'similarity'].median() for c in pairs.columns if c != 'similarity'],
        'different': [np.nan]
                     + [pairs.loc[~pairs[c], 'similarity'].median() for c in pairs.columns if c != 'similarity']})

    return table.round(3)


def draw(config, similarity, output_stem):
    matplotlib.rcParams['font.family'] = 'sans-serif'
    parasites = config['parasites']
    taxon_colors = config['parasite_groups']
    taxids = list(similarity.index)
    n = len(taxids)
    labels = [parasites[t]['label'] for t in taxids]

    fig, ax = plt.subplots(figsize=(7.0, 7.0))
    cmap = LinearSegmentedColormap.from_list('app', SCALE)
    # the diagonal compares nothing and is grey, not the white of a zero
    cmap.set_bad(DIAGONAL_COLOR)
    image = ax.imshow(np.ma.masked_invalid(similarity.to_numpy()), cmap=cmap, vmin=0,
                      vmax=np.nanmax(similarity.to_numpy()), interpolation='nearest')
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, fontsize=5.5, rotation=90, style='italic')
    ax.set_yticklabels(labels, fontsize=5.5, style='italic')
    ax.tick_params(length=0, pad=2)
    for side in ax.spines.values():
        side.set_visible(False)
    # a hairline grid between the cells, as the app's gaps
    ax.set_xticks(np.arange(-0.5, n), minor=True)
    ax.set_yticks(np.arange(-0.5, n), minor=True)
    ax.grid(which='minor', color='white', linewidth=0.4)
    ax.tick_params(which='minor', length=0)

    # the group and niche strips beside the rows and above the columns, a cell clear of
    # the matrix, so both axes visibly carry the same list in the same order
    strip = 0.8
    for i, t in enumerate(taxids):
        p = parasites[t]
        for offset, color in ((0, taxon_colors[p['group']]), (1, NICHE_COLORS[p['niche']])):
            start = -1.5 - strip - offset * strip * 1.1
            ax.add_patch(Rectangle((start, i - 0.5), strip * 0.9, 1, color=color,
                                   linewidth=0, clip_on=False))
            ax.add_patch(Rectangle((i - 0.5, start), 1, strip * 0.9, color=color,
                                   linewidth=0, clip_on=False))
    pad = 2 + (1 + 2.2 * strip) * (7.0 * 0.62 / n) * 72
    ax.tick_params(axis='x', pad=pad, top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.tick_params(axis='y', pad=pad)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)

    bar = fig.colorbar(image, ax=ax, fraction=0.03, pad=0.02)
    bar.ax.tick_params(labelsize=6, length=2)
    bar.set_label('Jaccard similarity of the host proteins reached', fontsize=6.5)
    bar.outline.set_visible(False)

    handles = [Patch(color=color, label=group) for group, color in taxon_colors.items()
               if any(parasites[t]['group'] == group for t in taxids)]
    handles += [Patch(color=color, label=niche) for niche, color in NICHE_COLORS.items()]
    fig.legend(handles=handles, loc='lower center', ncol=5, fontsize=6, frameon=False,
               handlelength=0.9, bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(left=0.26, right=0.9, top=0.78, bottom=0.07)
    for extension in ('pdf', 'svg'):
        fig.savefig(f'{output_stem}.{extension}')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--host', type=int, default=9606, help='host taxid (default: human)')
    parser.add_argument('--order-by', choices=['group', 'niche'], default='niche',
                        help='the annotation the parasites are ordered by first')
    parser.add_argument('--table', default='paper/tables/parasite_similarity.csv')
    parser.add_argument('--medians', default='paper/tables/parasite_similarity_medians.csv')
    parser.add_argument('--figure', default='paper/figures/parasite_similarity',
                        help='output path without extension')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    predictions = pd.read_parquet(os.path.join(args.data_dir, 'predictions.parquet'))
    similarity = similarity_matrix(interactors(config, predictions, args.host, args.order_by))
    named = similarity.rename(index=lambda t: config['parasites'][t]['label'],
                              columns=lambda t: config['parasites'][t]['label'])
    named.round(3).to_csv(args.table)
    table = medians(config, similarity)
    table.to_csv(args.medians, index=False)
    print(table.to_string(index=False))
    print(f"Wrote {args.table} and {args.medians}")
    draw(config, similarity, args.figure)
    print(f"Wrote {args.figure}.pdf and .svg")
