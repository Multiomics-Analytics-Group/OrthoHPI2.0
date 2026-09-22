'''
Which host protein families the parasites of one host converge on. A family is the
orthology group of the host protein, and a family is reached by a parasite when any of
its proteins is predicted to interact with any protein of the family. Beside the count
of parasites, the two figures that say whether the convergence is real: how many
distinct parasite orthology groups reach the family (one means a single conserved
parasite protein transferred many times), and how many helminths and protozoa do.
Writes paper/tables/shared_families.csv (every family reached by two or more parasites)
and paper/figures/shared_families.pdf/.svg (the families reached by at least
--min-parasites, as a dot matrix over the parasites).
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
from scripts import figure_style

HELMINTHS = {'Nematoda', 'Trematoda', 'Cestoda'}
NICHE_COLORS = {'extracellular': '#c7c7c7', 'intracellular': '#3d3d3d'}
DOT_COLOR = '#2a62ab'
# most gene symbols in a family label before the rest are left as an ellipsis
SYMBOLS_IN_LABEL = 2
LABEL_CHARS = 16


def family_label(names):
    '''The symbols of the family's proteins that were reached, the way the app names it.'''
    names = sorted(set(str(n).upper() for n in names))
    named = []
    for name in names[:SYMBOLS_IN_LABEL]:
        if named and len(', '.join(named + [name])) > LABEL_CHARS:
            break
        named.append(name)
    label = ', '.join(named)

    return f'{label}…' if len(named) < len(names) else label


def count_families(config, predictions, host):
    '''One row per host family reached by two or more parasites of the host.'''
    parasites = config['parasites']
    edges = predictions[predictions['taxid2'] == str(host)].copy()
    edges['taxon'] = edges['taxid1'].map(lambda t: parasites[int(t)]['group'])
    edges['clade'] = edges['taxon'].map(lambda g: 'helminth' if g in HELMINTHS else 'protozoan')

    families = edges.groupby('group2').agg(
        parasites=('taxid1', 'nunique'),
        helminths=('taxid1', lambda t: t[edges.loc[t.index, 'clade'] == 'helminth'].nunique()),
        protozoa=('taxid1', lambda t: t[edges.loc[t.index, 'clade'] == 'protozoan'].nunique()),
        parasite_groups=('group1', 'nunique'),
        host_proteins=('target', 'nunique'),
        interactions=('target', 'size'),
        symbols=('target_name', lambda n: ', '.join(sorted(set(str(x).upper() for x in n)))),
        label=('target_name', family_label))
    families = families[families['parasites'] > 1].sort_values(
        ['parasites', 'parasite_groups', 'interactions'], ascending=False, kind='stable')
    families.index.name = 'family'

    return families.reset_index(), edges


def draw(config, families, edges, output_stem):
    figure_style.apply()
    parasites = config['parasites']
    taxon_order = list(config['parasite_groups'])
    taxon_colors = config['parasite_groups']
    columns = sorted(set(edges['taxid1']),
                     key=lambda t: (taxon_order.index(parasites[int(t)]['group']),
                                    parasites[int(t)]['label']))
    col_of = {t: i for i, t in enumerate(columns)}
    row_of = {f: i for i, f in enumerate(families['family'])}

    # a dot per (family, parasite), its area the parasite proteins reaching the family
    dots = (edges[edges['group2'].isin(row_of)]
            .groupby(['group2', 'taxid1'])['source'].nunique().reset_index())
    n_rows, n_cols = len(row_of), len(columns)
    fig, ax = plt.subplots(figsize=(figure_style.WIDTH, figure_style.ROW * n_rows + 2.6))
    ax.scatter([col_of[t] for t in dots['taxid1']], [row_of[f] for f in dots['group2']],
               s=6 + 3 * dots['source'].clip(upper=20), color=DOT_COLOR, linewidths=0, zorder=3)
    ax.set_xlim(-0.6, n_cols - 0.4)
    ax.set_ylim(n_rows - 0.4, -0.6)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(families['label'])
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([figure_style.short_name(parasites[int(t)]['label']) for t in columns],
                       rotation=90, style='italic')
    ax.tick_params(length=0, pad=2)
    ax.grid(color='#eeeeee', linewidth=0.5)
    ax.set_axisbelow(True)
    for side in ax.spines.values():
        side.set_visible(False)
    # the parasite names go above the strips, the strips above the dots
    ax.xaxis.tick_top()

    # two strips over the columns: taxonomic group, then niche
    strip = 0.55
    for i, t in enumerate(columns):
        p = parasites[int(t)]
        ax.add_patch(Rectangle((i - 0.5, -0.6 - strip), 1, strip * 0.8,
                               color=taxon_colors[p['group']], linewidth=0, clip_on=False))
        ax.add_patch(Rectangle((i - 0.5, -0.6 - 2 * strip), 1, strip * 0.8,
                               color=NICHE_COLORS[p['niche']], linewidth=0, clip_on=False))
    ax.tick_params(axis='x', pad=2 + 2 * strip * figure_style.ROW * 72)

    # the margin columns: parasites, parasite groups, helminths / protozoa
    margin_x = n_cols + 0.6
    headers = [('parasites', 'parasites'), ('parasite_groups', 'parasite groups'),
               (None, 'helminths / protozoa')]
    for j, (column, header) in enumerate(headers):
        x = margin_x + 2.5 * j
        # the headers stand like the parasite names, from the top of the strips up
        ax.text(x, -0.6 - 2 * strip - 0.3, header, ha='center', va='bottom',
                rotation=90, color='#555555', clip_on=False)
        for row in families.itertuples():
            y = row_of[row.family]
            text = (f'{row.helminths} / {row.protozoa}' if column is None
                    else f'{getattr(row, column)}')
            ax.text(x, y, text, ha='center', va='center', color='#555555', clip_on=False)

    groups = [Patch(color=color, label=group) for group, color in taxon_colors.items()
              if any(parasites[int(t)]['group'] == group for t in columns)]
    niches = [Patch(color=color, label=niche) for niche, color in NICHE_COLORS.items()]
    legends = figure_style.stack_legends(fig, left=0.18, blocks=[('Taxonomic group', groups), ('Niche', niches)])
    # room above the dots for the strips and the parasite names standing on them
    fig.subplots_adjust(left=0.18, right=0.875, top=1 - 1.45 / fig.get_figheight(),
                        bottom=(legends + 0.1) / fig.get_figheight())
    for extension in ('pdf', 'svg'):
        fig.savefig(f'{output_stem}.{extension}')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--host', type=int, default=9606, help='host taxid (default: human)')
    parser.add_argument('--min-parasites', type=int, default=12,
                        help='families reached by fewer parasites stay out of the figure')
    parser.add_argument('--table', default='paper/tables/shared_families.csv')
    parser.add_argument('--figure', default='paper/figures/shared_families',
                        help='output path without extension')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    predictions = pd.read_parquet(os.path.join(args.data_dir, 'predictions.parquet'))
    families, edges = count_families(config, predictions, args.host)
    families.to_csv(args.table, index=False)
    print(f"Wrote {args.table}: {len(families)} families reached by two or more parasites, "
          f"{(families['parasites'] >= args.min_parasites).sum()} by {args.min_parasites} or more")
    draw(config, families[families['parasites'] >= args.min_parasites], edges, args.figure)
    print(f"Wrote {args.figure}.pdf and .svg")
