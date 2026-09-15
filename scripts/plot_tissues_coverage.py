'''
Bar plot per host of how many STRING proteins each config tissue has expression evidence
for: TISSUES (>= the host's cutoff) or the single-cell atlases: HPA for human (nTPM >
cutoff in some cell type) and Tabula Muris Senis / the pig atlas for mouse and pig, whose
values are mean log-normalised expression, so a gene counts where that is > 0.

Usage: .venv/bin/python scripts/plot_tissues_coverage.py [--source tissues|atlas] [--ntpm 1.0]
           [--output FILE]
'''
import argparse
import math
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import utils
from pipeline import cell_type_annotations, main as pipeline_main
import compare_tissue_filters
import venn_tissues_vs_hpa as venn

# a shared x-extent so the bars have the same width on every host
MAX_TISSUES = 22
BAR_COLOR = {'tissues': '#2a78d6', 'atlas': '#eb6834'}
INK = '#0b0b0b'
INK_SECONDARY = '#52514e'


def atlas_tissue_proteins(config_file, data_dir, taxid, pool, ntpm_cutoff):
    '''{config tissue label: pool proteins above the cutoff in some cell type}, per atlas.'''
    if taxid == venn.HUMAN:
        data = cell_type_annotations.read_hpa(config_file=config_file)
        data = cell_type_annotations.map_hpa_data(config_file=config_file, hpa_data=data)
    elif str(taxid) == cell_type_annotations.MOUSE_TAXID:
        data = cell_type_annotations.read_mouse_atlas(data_dir, pool)
    elif str(taxid) == cell_type_annotations.PIG_TAXID:
        data = cell_type_annotations.read_pig_atlas(data_dir, pool)
    else:
        return {}
    cutoff = ntpm_cutoff if taxid == venn.HUMAN else 0.0
    data = data[data['Gene'].isin(pool) & (data['nTPM'] > cutoff)]
    return data.groupby('Tissue')['Gene'].apply(set).to_dict()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--source', choices=['tissues', 'atlas'], default='tissues')
    parser.add_argument('--ntpm', type=float, default=venn.NTPM_CUTOFF, help='atlas nTPM cutoff')
    parser.add_argument('--output', help='PNG to write (default: snapshots/<source>_coverage.png)')
    args = parser.parse_args()
    output = args.output or os.path.join('snapshots', f'{args.source}_coverage.png')

    hosts = utils.read_config(filepath=args.config, field='hosts')
    pools = compare_tissue_filters.host_pools(args.config)
    rows = math.ceil(len(hosts) / 2)
    fig, axes = plt.subplots(rows, 2, figsize=(16, 4.2 * rows))
    for ax in axes.flat[len(hosts):]:
        ax.axis('off')
    for ax, (taxid, host) in zip(axes.flat, hosts.items()):
        pool = pools[taxid]
        if args.source == 'tissues':
            by_tissue = venn.raw_tissues_tissue_proteins(args.config, pool, taxid)
            ylabel = f"proteins, TISSUES score ≥ {host.get('tissue_cutoff', pipeline_main.TISSUE_CUTOFF):g}"
        else:
            by_tissue = atlas_tissue_proteins(args.config, args.data_dir, taxid, pool, args.ntpm)
            ylabel = (f'proteins, nTPM > {args.ntpm:g} in a cell type' if taxid == venn.HUMAN
                      else 'proteins detected in a cell type')
        counts = sorted(((len(p), t) for t, p in by_tissue.items()), reverse=True)
        if not counts:
            ax.text(0.5, 0.5, 'no single-cell atlas', ha='center', va='center', transform=ax.transAxes,
                    fontsize=12, color=INK_SECONDARY)
            ax.set_title(f"{host['label']}", fontsize=11, color=INK, loc='left')
            ax.axis('off')
            continue
        top = max(n for n, _ in counts)

        ax.bar([t for _, t in counts], [n for n, _ in counts], color=BAR_COLOR[args.source], width=0.72)
        for x, (n, _) in enumerate(counts):
            ax.text(x, n + 0.01 * top, f'{n:,}', ha='center', fontsize=7, color=INK)
        ax.set_xlim(-0.6, MAX_TISSUES - 0.4)
        ax.set_ylim(0, top * 1.1)
        ax.set_ylabel(ylabel, color=INK_SECONDARY, fontsize=9)
        ax.tick_params(axis='x', length=0, labelsize=8, rotation=45)
        plt.setp(ax.get_xticklabels(), ha='right', rotation_mode='anchor')
        ax.tick_params(axis='y', colors=INK_SECONDARY, labelsize=8)
        for spine in ('top', 'right', 'bottom'):
            ax.spines[spine].set_visible(False)
        ax.spines['left'].set_color('#c9c8c2')
        ax.set_title(f"{host['label']}: {len(counts)} tissues annotated, of {len(pool):,} STRING proteins",
                     fontsize=11, color=INK, loc='left')
    fig.tight_layout(h_pad=2.5, w_pad=2)

    os.makedirs(os.path.dirname(output) or '.', exist_ok=True)
    fig.savefig(output, dpi=200, bbox_inches='tight', facecolor='white')
    print('wrote', output)


if __name__ == '__main__':
    main()
