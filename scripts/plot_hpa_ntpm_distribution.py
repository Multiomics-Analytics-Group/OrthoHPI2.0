'''
Ridgeline of the HPA single-cell nTPM distribution, one ridge per tissue (--by tissue)
or per cell type (--by cell-type), over the tissues and genes OrthoHPI maps.

nTPM spans five orders of magnitude, so the ridges are histograms of log10(nTPM) over the
detected values; the undetected (nTPM = 0) share is printed in the left gutter and counts
towards the "% kept" columns on the right, which are the share of all (gene, tissue, cell
type) combinations a cutoff would keep.

--per-gene-max reduces each (gene, tissue) to its highest cell type first, which is what the
pipeline's HPA filter actually tests.

Usage: .venv/bin/python scripts/plot_hpa_ntpm_distribution.py [--by tissue|cell-type]
           [--tissue TISSUE] [--per-gene-max] [--cutoffs 1 20] [--columns 3] [--output FILE]
'''
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib.pyplot as plt
import utils
from pipeline import cell_type_annotations

# dataviz palette: one measure, so one hue
FILL_COLOR = '#2a78d6'
CUTOFF_COLOR = '#eb6834'
INK = '#0b0b0b'
INK_SECONDARY = '#52514e'
INK_MUTED = '#8a8880'
GRID = '#dcdbd5'


def read_hpa_with_zeros(config_file):
    '''The HPA single-cell table mapped to config tissues and STRING proteins, zeros kept.'''
    urls = utils.read_config(filepath=config_file, field='urls')
    filename = utils.download_file(url=urls['hpa_single_cell_tissue_url'], data_dir='data/downloads')
    data = pd.read_csv(utils.read_zipped_file(filepath=filename), sep='\t', header=0)
    data = data.sort_values(by='nTPM', ascending=False)
    data = data.drop_duplicates(['Gene', 'Tissue', 'Cell type'], keep='first')
    data = cell_type_annotations.map_hpa_data(config_file=config_file, hpa_data=data)
    return data.dropna(subset=['Gene'])


def ridge_rows(data, group, cutoffs):
    '''One record per group: its detected log10 values, its zero share and its kept shares.'''
    rows = []
    for name, values in data.groupby(group)['nTPM']:
        detected = values[values > 0]
        rows.append({
            'name': name,
            'n': len(values),
            'log10': np.log10(detected.to_numpy()),
            'zero': 1 - len(detected) / len(values),
            'median': detected.median() if len(detected) else 0.0,
            'kept': [(values > c).mean() for c in cutoffs],
        })
    return sorted(rows, key=lambda row: row['median'])


def draw_panel(ax, rows, cutoffs, lo, hi, edges, centres, kernel, overlap, ymax, gutters=True,
               label_cutoffs=True):
    '''One column of ridges, with the cutoff rules and, unless suppressed, the share columns.'''
    for cutoff in cutoffs:
        ax.axvline(np.log10(cutoff), color=CUTOFF_COLOR, linewidth=1.2, linestyle=(0, (4, 3)),
                   zorder=1)
    for x in range(int(lo), int(hi) + 1):
        ax.axvline(x, color=GRID, linewidth=0.6, zorder=0)

    for i, row in enumerate(rows):
        density, _ = np.histogram(np.clip(row['log10'], lo, hi), bins=edges, density=True)
        density = np.convolve(density, kernel, mode='same')
        peak = density.max() or 1.0
        y = i + overlap * density / peak
        ax.fill_between(centres, i, y, color=FILL_COLOR, alpha=0.55, linewidth=0, zorder=2 + i)
        ax.plot(centres, y, color=FILL_COLOR, linewidth=1.4, zorder=2 + i)
        ax.plot([lo, hi], [i, i], color='white', linewidth=2, zorder=1.9 + i)

    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([row['name'] for row in rows], fontsize=9, color=INK)
    ax.set_ylim(-0.6, ymax + overlap)
    ax.set_xlim(lo, hi)
    ax.set_xticks(range(int(lo), int(hi) + 1))
    ax.set_xticklabels([f'{10.0 ** x:g}' for x in range(int(lo), int(hi) + 1)], color=INK_SECONDARY)
    ax.set_xlabel('nTPM (log scale)', color=INK_SECONDARY)
    ax.tick_params(axis='both', length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # the cutoff rules label themselves at the top of the plot, where a panel is wide
    # enough; in a narrow one the labels would run into each other and the subtitle names them
    for cutoff in cutoffs if label_cutoffs else ():
        ax.annotate(f'nTPM > {cutoff:g}', xy=(np.log10(cutoff), ymax + overlap),
                    xytext=(0, 4), textcoords='offset points', ha='center', va='bottom',
                    fontsize=9, color=CUTOFF_COLOR, fontweight='bold', annotation_clip=False,
                    zorder=100)

    if not gutters:
        return
    # the undetected share, then the share each cutoff keeps
    columns = [('zero', [row['zero'] for row in rows], INK_MUTED)]
    columns += [(f'> {cutoff:g}', [row['kept'][j] for row in rows], INK)
                for j, cutoff in enumerate(cutoffs)]
    for j, (header, values, color) in enumerate(columns):
        x = 18 + 46 * j
        ax.annotate(header, xy=(1, 1), xytext=(x, 6), textcoords='offset points',
                    xycoords='axes fraction', ha='center', va='bottom', fontsize=8.5,
                    color=INK_MUTED)
        for i, value in enumerate(values):
            ax.annotate(f'{value:.0%}', xy=(1, i), xytext=(x, 0), textcoords='offset points',
                        xycoords=('axes fraction', 'data'), ha='center', va='center',
                        fontsize=8.5, color=color, annotation_clip=False)


def ridge_axis(rows):
    '''The shared log10 range, bins and smoothing kernel every panel is drawn on.'''
    lo = -1
    # the long thin tail past the 99.9th percentile would leave most of the axis empty
    hi = np.ceil(np.percentile(np.concatenate([row['log10'] for row in rows]), 99.9))
    edges = np.linspace(lo, hi, 160)
    centres = (edges[:-1] + edges[1:]) / 2
    # nTPM is reported to one decimal, so the raw histogram combs below ~1; smooth it
    kernel = np.exp(-0.5 * (np.arange(-12, 13) / 3.0) ** 2)
    kernel /= kernel.sum()
    return lo, hi, edges, centres, kernel


def draw_facets(groups, cutoffs, title, subtitle, output, columns=4):
    '''A grid of panels, one per tissue, each holding that tissue's cell types as ridges.'''
    lo, hi, edges, centres, kernel = ridge_axis([row for _, rows in groups for row in rows])
    overlap = 0.92

    # biggest panels first, so each grid row is only as tall as it has to be
    groups = sorted(groups, key=lambda g: -len(g[1]))
    grid = [groups[i:i + columns] for i in range(0, len(groups), columns)]
    heights = [max(len(rows) for _, rows in band) for band in grid]

    fig = plt.figure(figsize=(5.0 * columns + 1.0, 0.30 * sum(heights) + 1.3 * len(grid) + 1.6))
    spec = fig.add_gridspec(len(grid), columns, height_ratios=heights, hspace=0.55, wspace=0.75)
    for r, band in enumerate(grid):
        for c, (name, rows) in enumerate(band):
            ax = fig.add_subplot(spec[r, c])
            draw_panel(ax, rows, cutoffs, lo, hi, edges, centres, kernel, overlap,
                       heights[r], gutters=False, label_cutoffs=False)
            ax.set_title(name, fontsize=10.5, color=INK, fontweight='bold', loc='left', x=0, pad=16)

    fig.suptitle(f'{title}\n{subtitle}', fontsize=11, color=INK, x=0.02, ha='left')
    fig.savefig(output, dpi=200, bbox_inches='tight', facecolor='white')
    return output


def draw(rows, cutoffs, title, subtitle, output, columns=1):
    '''Ridgeline over log10(nTPM). Several columns split the rows over side-by-side panels.'''
    lo, hi, edges, centres, kernel = ridge_axis(rows)
    # side by side there is room to keep every ridge inside its own band; stacked in one
    # column they have to bleed over their neighbours to stay readable
    overlap = 0.92 if columns > 1 else (1.9 if len(rows) <= 25 else 1.5)

    # highest median first, filling each column from the top before starting the next
    panels = [list(reversed(list(chunk)))
              for chunk in np.array_split(list(reversed(rows)), columns)]
    tallest = max(len(panel) for panel in panels)

    if columns == 1:
        fig, ax = plt.subplots(figsize=(11, 0.34 * len(rows) + 2.6))
        axes = [ax]
    else:
        # each panel needs room for its own labels plus the share columns of the one before
        fig, axes = plt.subplots(1, columns, squeeze=False,
                                 figsize=(7.8 * columns + 1.4, 0.34 * tallest + 2.6))
        axes = list(axes[0])
        fig.subplots_adjust(wspace=1.1)

    for ax, panel in zip(axes, panels):
        draw_panel(ax, panel, cutoffs, lo, hi, edges, centres, kernel, overlap, tallest)

    if columns == 1:
        axes[0].set_title(f'{title}\n{subtitle}', fontsize=11, color=INK, pad=38, loc='left', x=0)
    else:
        fig.suptitle(f'{title}\n{subtitle}', fontsize=11, color=INK, x=0.02, ha='left')

    fig.savefig(output, dpi=200, bbox_inches='tight', facecolor='white')
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--by', choices=['tissue', 'cell-type', 'tissue-cell-type'],
                        default='tissue',
                        help='one ridge per config tissue, per HPA cell type, or, for '
                             'tissue-cell-type, a grid of one panel per tissue holding its own '
                             'cell types (default: tissue)')
    parser.add_argument('--tissue', help='with --by cell-type: only this tissue\'s cell types')
    parser.add_argument('--per-gene-max', action='store_true',
                        help='reduce each (gene, tissue) to its highest cell type, as the filter does')
    parser.add_argument('--columns', type=int, default=1,
                        help='split the rows over this many side-by-side panels, which leaves '
                             'room for the ridges not to overlap (default: 1)')
    parser.add_argument('--cutoffs', type=float, nargs='*', default=[1.0, 20.0],
                        help='nTPM cutoffs to mark (default: 1 20)')
    parser.add_argument('--output', help='PNG to write (default: snapshots/hpa_ntpm_<by>.png)')
    args = parser.parse_args()
    if args.per_gene_max and args.by != 'tissue':
        parser.error('--per-gene-max collapses the cell types, so it only applies to --by tissue')

    data = read_hpa_with_zeros(args.config)
    if args.tissue:
        if args.by != 'cell-type':
            parser.error('--tissue only applies to --by cell-type')
        data = data[data['Tissue'] == args.tissue.lower()]
        if data.empty:
            parser.error(f'no HPA data for tissue {args.tissue!r}')
    if args.per_gene_max:
        data = data.groupby(['Gene', 'Tissue'], as_index=False)['nTPM'].max()

    unit = 'gene, tissue' if args.per_gene_max else 'gene, tissue, cell type'
    output = args.output or os.path.join('snapshots', f"hpa_ntpm_{args.by.replace('-', '')}"
                                                      f"{'_pergenemax' if args.per_gene_max else ''}.png")
    os.makedirs(os.path.dirname(output) or '.', exist_ok=True)

    def report(rows, prefix=()):
        return [{**dict(zip(('tissue',), prefix)), 'name': row['name'], 'n': row['n'],
                 'zero %': round(100 * row['zero'], 1), 'median nTPM': row['median'],
                 **{f'% > {c:g}': round(100 * k, 1) for c, k in zip(args.cutoffs, row['kept'])}}
                for row in reversed(rows)]

    if args.by == 'tissue-cell-type':
        groups = [(tissue, ridge_rows(part, 'Cell type', args.cutoffs))
                  for tissue, part in data.groupby('Tissue')]
        title = 'How HPA single-cell nTPM is distributed across the cell types of each tissue'
        subtitle = (f'{len(data):,} ({unit}) combinations over '
                    f'{sum(len(rows) for _, rows in groups)} (tissue, cell type) pairs in '
                    f'{len(groups)} tissues; ridges are the detected values, cell types '
                    'ordered by median within each tissue; dashed rules mark nTPM = '
                    + ' and '.join(f'{c:g}' for c in args.cutoffs))
        print(pd.DataFrame([r for tissue, rows in groups
                            for r in report(rows, (tissue,))]).to_string(index=False))
        columns = args.columns if args.columns > 1 else 4
        print('\nwrote', draw_facets(groups, args.cutoffs, title, subtitle, output, columns))
        return

    group = 'Tissue' if args.by == 'tissue' else 'Cell type'
    rows = ridge_rows(data, group, args.cutoffs)
    noun = 'tissues' if args.by == 'tissue' else 'cell types'
    title = f'How HPA single-cell nTPM is distributed across {noun}'
    subtitle = (f'{len(data):,} ({unit}) combinations over {len(rows)} {noun}'
                + (f" in {args.tissue.lower()}" if args.tissue else '')
                + ('; each gene at its highest cell type' if args.per_gene_max else '')
                + '; ridges are the detected values, rows ordered by median')
    print(pd.DataFrame(report(rows)).to_string(index=False))
    print('\nwrote', draw(rows, args.cutoffs, title, subtitle, output, args.columns))


if __name__ == '__main__':
    main()
