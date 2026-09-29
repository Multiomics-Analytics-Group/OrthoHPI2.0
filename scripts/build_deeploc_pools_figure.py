'''
The DeepLoc 2 probabilities of every protein a localization filter decided on, kept or
not, so the cut-off is seen against the whole distribution it cuts. Host proteins: one
panel per class, a violin per host over the proteins expressed in a tissue infected by
one of its parasites whose niche reads that class (the tissue filter comes first, so
these are the proteins the threshold decides on). Parasite proteins: the whole proteome
of each parasite, on extracellular, and on cell membrane for the unicellular parasites,
the secretome filter reading that class for them alone. The part of a violin over the
cut-off is shaded, and the share of the pool it keeps is printed beside it. Writes
paper/tables/deeploc_pools.csv (one row per violin, with the proteins within 0.05 either
side of the cut-off) and paper/figures/deeploc_pools_host.svg/.png and
deeploc_pools_parasite.svg/.png.
'''
import argparse
import copy
import glob
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils
from pipeline import filters, main
from scripts import figure_style
from scripts.build_deeploc_figure import CLASS_COLORS, NICHE_COLORS, SURFACE

# how far either side of a cut-off counts as near it in the table
MARGIN = 0.05


def read_deeploc(data_dir, taxid):
    '''The DeepLoc results of one species, indexed by STRING id.'''
    matches = sorted(glob.glob(os.path.join(data_dir, main.DEEPLOC_ACCURATE_DIR, str(taxid),
                                            'results_*.csv')))

    return pd.read_csv(matches[-1]).set_index('Protein_ID')


def host_pools(config_file, config, data_dir):
    '''{(host, class): probabilities of the host proteins the threshold of that class
    decided on}; hosts in config order.'''
    proteins = main.get_proteins(config_file)
    tissues = filters.apply_tissue_filter(config_file=config_file,
                                          valid_proteins=copy.deepcopy(proteins),
                                          cutoff=main.TISSUE_CUTOFF)
    infected = filters.parasite_tissue_proteins(config_file, tissues, proteins)

    pools = {}
    for host_taxid, host in config['hosts'].items():
        name = host['label'].split('(')[1].rstrip(')')
        deeploc = read_deeploc(data_dir, host_taxid)
        prefix = f'{host_taxid}.'
        for surface_class in main.DEEPLOC_CUTOFFS:
            pool = set()
            for taxid, parasite in config['parasites'].items():
                if host_taxid in parasite['hosts'] and \
                        surface_class in main.DEEPLOC_NICHE_CLASSES[parasite['niche']]:
                    pool |= {p for p in infected[int(taxid)] if p.startswith(prefix)}
            pools[(name, surface_class)] = deeploc.loc[deeploc.index.isin(pool), surface_class]

    return pools


def parasite_pools(config, data_dir):
    '''{(parasite, class): probabilities of its whole proteome}, on extracellular for
    every parasite and cell membrane for the unicellular ones; parasites by group, niche
    then name.'''
    groups = list(config['parasite_groups'])
    ordered = sorted(config['parasites'].items(),
                     key=lambda item: (groups.index(item[1]['group']),
                                       list(NICHE_COLORS).index(item[1]['niche']),
                                       item[1]['label']))
    pools = {}
    for taxid, parasite in ordered:
        deeploc = read_deeploc(data_dir, taxid)
        for surface_class in SURFACE:
            if surface_class == 'Cell membrane' and parasite.get('multicellular', False):
                continue
            pools[(parasite['label'], surface_class)] = deeploc[surface_class]

    return pools


def summarise(host, parasite):
    '''One row per violin: the pool, what the cut-off keeps of it, and how many proteins
    sit within MARGIN under and over it.'''
    rows = []
    for side, pools in (('host', host), ('parasite', parasite)):
        for (species, surface_class), values in pools.items():
            cutoff = main.DEEPLOC_CUTOFFS[surface_class]
            kept = int((values > cutoff).sum())
            rows.append({'side': side, 'localization': surface_class, 'species': species,
                         'pool': len(values), 'kept': kept,
                         'kept_share': kept / len(values) if len(values) else float('nan'),
                         'under_cutoff': int(((values > cutoff - MARGIN) & (values <= cutoff)).sum()),
                         'over_cutoff': int(((values > cutoff) & (values <= cutoff + MARGIN)).sum())})

    return pd.DataFrame(rows)


def draw_violin(ax, values, position, color, cutoff):
    '''A violin over `values`, pale under the cut-off and full over it, with the share
    it keeps at the right.'''
    for alpha, clip in ((0.3, None), (0.9, cutoff)):
        body = ax.violinplot(values, positions=[position], vert=False, widths=0.85,
                             showextrema=False, points=200)['bodies'][0]
        body.set_facecolor(color)
        body.set_edgecolor('none')
        body.set_alpha(alpha)
        if clip is not None:
            body.set_clip_path(Rectangle((clip, position - 1), 1, 2, transform=ax.transData))
    ax.text(1.02, position, f'{(values > cutoff).mean():.0%}', va='center', color='#444444',
            transform=ax.get_yaxis_transform())


def style_axis(ax, cutoff, n_rows):
    ax.axvline(cutoff, color='#333333', linestyle=':', linewidth=1, zorder=3)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_ylim(n_rows - 0.4, -0.6)
    ax.tick_params(axis='y', length=0)
    ax.grid(axis='x', color='#e6e6e6', linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_color('#999999')


def draw_hosts(pools, table, output_stem):
    '''Four panels, one per class, a row per host.'''
    figure_style.apply()
    hosts = list(dict.fromkeys(host for host, _ in pools))
    fig, axes = plt.subplots(2, 2, figsize=(figure_style.WIDTH, 6.4), sharey=True)
    for letter, ax, (surface_class, cutoff) in zip('ABCD', axes.flat, main.DEEPLOC_CUTOFFS.items()):
        for i, host in enumerate(hosts):
            values = pools[(host, surface_class)]
            if len(values):
                draw_violin(ax, values, i, CLASS_COLORS[surface_class], cutoff)
            else:
                ax.text(0.01, i, 'not read', va='center', color='#888888')
        style_axis(ax, cutoff, len(hosts))
        ax.set_yticks(range(len(hosts)))
        ax.set_yticklabels(hosts)
        ax.set_title(f'({letter}) {surface_class}', fontweight='bold', loc='left', pad=3)
        ax.set_xlabel(f'P({surface_class.lower()})')
    fig.tight_layout(h_pad=1.2, w_pad=3)
    for extension in ('svg', 'png'):
        fig.savefig(f'{output_stem}.{extension}', dpi=300)
    plt.close(fig)


def draw_parasites(config, pools, output_stem):
    '''Two panels side by side, every parasite on extracellular and the unicellular ones
    on cell membrane, on a shared row grid.'''
    figure_style.apply()
    taxon_colors = config['parasite_groups']
    by_label = {p['label']: p for p in config['parasites'].values()}
    rows = {c: [p for p, pc in pools if pc == c] for c in SURFACE}
    n_rows = max(len(r) for r in rows.values())
    foot = 1.55
    row = min(figure_style.ROW, (figure_style.MAX_HEIGHT - foot - 0.35) / n_rows)
    fig = plt.figure(figsize=(figure_style.WIDTH, row * n_rows + foot + 0.35))
    grid = fig.add_gridspec(n_rows, 2, wspace=1.6,
                            top=1 - 0.35 / fig.get_figheight(), bottom=foot / fig.get_figheight(),
                            left=0.19, right=0.93)
    for column, (letter, surface_class) in enumerate(zip('AB', SURFACE)):
        names = rows[surface_class]
        ax = fig.add_subplot(grid[:len(names), column])
        cutoff = main.DEEPLOC_CUTOFFS[surface_class]
        for i, name in enumerate(names):
            parasite = by_label[name]
            draw_violin(ax, pools[(name, surface_class)], i, taxon_colors[parasite['group']], cutoff)
            # the group and niche strips beside the names
            for offset, color in ((0, taxon_colors[parasite['group']]), (1, NICHE_COLORS[parasite['niche']])):
                ax.add_patch(Rectangle((-0.8 - 0.035 * offset, i - 0.4), 0.022, 0.8, color=color,
                                       linewidth=0, clip_on=False, transform=ax.get_yaxis_transform()))
        style_axis(ax, cutoff, len(names))
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels([figure_style.short_name(p) for p in names], style='italic')
        ax.set_title(f'({letter}) {surface_class}', fontweight='bold', loc='left', pad=3)
        ax.set_xlabel(f'P({surface_class.lower()})')

    figure_style.stack_legends(fig, left=0.19, blocks=[
        ('Taxonomic group', [Patch(color=c, label=g) for g, c in taxon_colors.items()]),
        ('Niche', [Patch(color=c, label=niche) for niche, c in NICHE_COLORS.items()])])
    for extension in ('svg', 'png'):
        fig.savefig(f'{output_stem}.{extension}', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--table', default='paper/tables/deeploc_pools.csv')
    parser.add_argument('--figures', default='paper/figures/deeploc_pools',
                        help='output path without the _host/_parasite suffix and extension')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    host = host_pools(args.config, config, args.data_dir)
    parasite = parasite_pools(config, args.data_dir)
    table = summarise(host, parasite)
    table.to_csv(args.table, index=False, float_format='%.3f')
    print(f"{len(table)} pools; wrote {args.table}")
    draw_hosts(host, table, f'{args.figures}_host')
    draw_parasites(config, parasite, f'{args.figures}_parasite')
    print(f"Wrote {args.figures}_host and {args.figures}_parasite .svg and .png")
