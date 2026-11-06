'''
Which host orthologous groups the parasites of one host share across niches. A group is
reached by a parasite when any of its proteins is predicted to interact with any protein
of the group, and is shared when at least --min-share of the parasites of each niche reach
it, so that the niche with more parasites does not decide the selection. Beside the counts
per niche, the figure that says whether the convergence is real: how many distinct parasite
orthologous groups reach the group (one means a single conserved parasite protein
transferred many times).
Writes paper/tables/shared_ogs.csv (every group reached by two or more parasites, with
how many helminths, protozoa, intracellular and extracellular parasites reach it) and
paper/figures/shared_ogs.svg/.png (the shared groups, as a dot matrix over the parasites).

Convergence within a genus is the other end of the same question, and too small to
draw: paper/tables/shared_ogs_genus.csv holds, for every genus with several species of
the host, the groups reached by all of its species and by no other parasite, with the
parasite groups the links came from, so that a convergence driven by one conserved
parasite protein can be told from one that is not.
'''
import argparse
import os
import sys
import textwrap

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils
from scripts import figure_style

HELMINTHS = {'Nematoda', 'Trematoda', 'Cestoda'}
NICHE_COLORS = {'intracellular': '#3d3d3d', 'extracellular': '#c7c7c7'}
DOT_COLOR = '#2a62ab'
# group labels wrap onto lines of at most this many characters, which the longest
# family name fits on
LABEL_CHARS = 40
# The family each group of the figure belongs to, so that a row can be found from the
# name the text uses rather than from its gene symbols. A group with no name here is
# labelled by its symbols alone.
GROUP_NAMES = {
    'KOG0994': 'Laminins',
    'KOG1226': 'Integrin β subunits',
    'KOG3540': 'Amyloid precursor proteins',
    'KOG2369': 'LCAT/PLA2G15 acyltransferases',
    'KOG0782': 'Diacylglycerol kinases',
    'KOG3058': 'Sphingomyelin synthases',
    'KOG4160': 'Lipid-transfer proteins',
    'KOG4338': 'Apolipoprotein B',
    'KOG3781': 'Dystroglycan',
    'KOG1450': 'Rho GTPase-activating proteins',
    'KOG3603': 'Phospholipase D',
    'KOG1480': 'Netrin receptors UNC5',
    'KOG1348': 'Legumain',
    'KOG4119': 'G protein γ subunits',
    'KOG1025': 'EGF receptor family',
    'KOG4292': 'CUB-domain proteins',
    'KOG1544': 'Tubulointerstitial nephritis antigens',
    'KOG3108': 'OB-fold proteins',
    'KOG4087': 'Secreted phospholipases A2',
    'KOG3509': 'Basement-membrane proteoglycans',
    'KOG3724': 'Uncharacterized C2ORF66',
    'KOG4221': 'Netrin receptors DCC/neogenin',
    'KOG1933': 'Sterol transporters NPC1/NPC1L1',
    'KOG2091': 'Chitinase-like CHID1',
    'KOG4257': 'Focal adhesion kinases',
    'KOG3065': 'SNAP25-family SNAREs',
    'KOG2497': 'Methyltransferase-like proteins',
    'KOG3511': 'Sortilins',
    'KOG0099': 'G protein α subunits',
    'KOG1935': 'Patched receptors',
    'KOG4004': 'SPARC/osteonectin',
    'KOG0264': 'Histone-binding RBBP4/RBBP7',
    'KOG1654': 'ATG8 family (LC3/GABARAP)',
    'KOG0366': 'Rab geranylgeranyltransferase β',
    'KOG0595': 'ULK kinases',
    'KOG0668': 'Casein kinase 2 α subunits',
    'KOG0603': 'Ribosomal S6 kinases',
    'KOG3092': 'Casein kinase 2 β subunit',
    'KOG1786': 'BEACH-domain proteins',
    'KOG2999': 'ELMO proteins',
    'KOG1091': 'Rab GTPase-activating TBC1D5',
    'KOG3874': 'Autophagy protein ATG13',
    'KOG4572': 'Autophagy protein RB1CC1/FIP200',
    'KOG1998': 'DOCK exchange factors',
    'KOG2171': 'Importins IPO4/IPO5/RANBP6',
    'KOG0285': 'Spliceosome protein PLRG1',
    'KOG0289': 'Spliceosome protein PRPF19',
    'KOG0282': 'Spliceosome protein CDC40',
    'KOG0643': 'Translation initiation factor EIF3I',
    'KOG0860': 'VAMP SNAREs',
    'KOG1299': 'Sec1-family VPS45',
    'KOG3251': 'Golgi SNARE GOSR2',
    'KOG1302': 'Sec1-family VPS33A/VPS33B',
    'KOG0834': 'Cyclins K, Q and T',
    'KOG2066': 'HOPS subunit VPS41',
    'KOG2063': 'HOPS/CORVET VPS39/VPS3',
    'KOG0811': 'Syntaxins STX7/STX12',
    'KOG0085': 'G protein αq GNA14/GNA15',
}
# the left margin, a figure fraction, that the wrapped labels take
LEFT = 0.27


def count_ogs(config, predictions, host):
    '''One row per host orthologous group reached by two or more parasites of the host.'''
    parasites = config['parasites']
    edges = predictions[predictions['taxid2'] == str(host)].copy()
    edges['taxon'] = edges['taxid1'].map(lambda t: parasites[int(t)]['group'])
    edges['clade'] = edges['taxon'].map(lambda g: 'helminth' if g in HELMINTHS else 'protozoan')
    edges['niche'] = edges['taxid1'].map(lambda t: parasites[int(t)]['niche'])

    ogs = edges.groupby('group2').agg(
        parasites=('taxid1', 'nunique'),
        helminths=('taxid1', lambda t: t[edges.loc[t.index, 'clade'] == 'helminth'].nunique()),
        protozoa=('taxid1', lambda t: t[edges.loc[t.index, 'clade'] == 'protozoan'].nunique()),
        intracellular=('taxid1', lambda t: t[edges.loc[t.index, 'niche'] == 'intracellular'].nunique()),
        extracellular=('taxid1', lambda t: t[edges.loc[t.index, 'niche'] == 'extracellular'].nunique()),
        parasite_groups=('group1', 'nunique'),
        host_proteins=('target', 'nunique'),
        interactions=('target', 'size'),
        symbols=('target_name', lambda n: ', '.join(sorted(set(str(x).upper() for x in n)))))
    ogs = ogs[ogs['parasites'] > 1].sort_values(
        ['parasites', 'parasite_groups', 'interactions'], ascending=False, kind='stable')
    ogs.index.name = 'orthologous_group'

    return ogs.reset_index(), edges


def count_genus_ogs(config, edges):
    '''
    One row per (genus, host orthologous group) reached by every species of the genus
    that has predictions, and by no parasite outside it, for the genera with several
    such species of the host.
    '''
    parasites = config['parasites']
    edges = edges.assign(species=edges['taxid1'].map(lambda t: parasites[int(t)]['label']))
    edges = edges.assign(genus=edges['species'].str.split().str[0])
    reached = edges.groupby('group2')['species'].apply(set)

    rows = []
    for genus, species in sorted(edges.groupby('genus')['species'].apply(set).items()):
        if len(species) < 2:
            continue
        exclusive = [group for group, found in reached.items() if found == species]
        for group in exclusive:
            links = edges[(edges['group2'] == group) & (edges['genus'] == genus)]
            scores = links['weight'].astype(float)
            rows.append({
                'genus': genus,
                'species': len(species),
                'orthologous_group': group,
                'family': GROUP_NAMES.get(group, ''),
                'host_proteins': links['target'].nunique(),
                'symbols': ', '.join(sorted(set(str(n).upper() for n in links['target_name']))),
                'parasite_groups': links['group1'].nunique(),
                'parasite_group_ids': ', '.join(sorted(set(links['group1']))),
                'parasite_proteins': links['source'].nunique(),
                'interactions': len(links),
                'min_score': round(scores.min(), 3),
                'max_score': round(scores.max(), 3)})

    columns = ['genus', 'species', 'orthologous_group', 'family', 'host_proteins', 'symbols',
               'parasite_groups', 'parasite_group_ids', 'parasite_proteins', 'interactions',
               'min_score', 'max_score']
    table = pd.DataFrame(rows, columns=columns)

    return table.sort_values(['genus', 'interactions'], ascending=[True, False], kind='stable')


def shared_ogs(ogs, edges, share):
    '''The groups reached by at least `share` of the parasites of each niche.'''
    sizes = edges.groupby('niche')['taxid1'].nunique()

    return ogs[(ogs['intracellular'] >= share * sizes['intracellular'])
               & (ogs['extracellular'] >= share * sizes['extracellular'])]


def row_label(group, symbols):
    '''
    The family name of the group, the name the text calls it by. The members reached
    are in the table, which holds too many symbols to label the rows with. A group with
    no name falls back to them.
    '''
    return textwrap.fill(GROUP_NAMES.get(group) or symbols, LABEL_CHARS)


def draw(config, ogs, edges, output_stem):
    figure_style.apply()
    parasites = config['parasites']
    taxon_order = list(config['parasite_groups'])
    taxon_colors = config['parasite_groups']
    columns = sorted(set(edges['taxid1']),
                     key=lambda t: (list(NICHE_COLORS).index(parasites[int(t)]['niche']),
                                    taxon_order.index(parasites[int(t)]['group']),
                                    parasites[int(t)]['label']))
    col_of = {t: i for i, t in enumerate(columns)}
    labels = [row_label(group, symbols)
              for group, symbols in zip(ogs['orthologous_group'], ogs['symbols'])]
    ys = figure_style.row_positions(labels, figure_style.DENSE_ROW)
    row_of = dict(zip(ogs['orthologous_group'], ys))

    # a dot per (group, parasite), its area the parasite proteins reaching the group
    dots = (edges[edges['group2'].isin(row_of)]
            .groupby(['group2', 'taxid1'])['source'].nunique().reset_index())
    n_cols = len(columns)
    fig, ax = plt.subplots(figsize=(figure_style.WIDTH, figure_style.DENSE_ROW * (ys[-1] + 1) + 2.35))
    size = lambda n: 6 + 3 * min(n, 20)
    ax.scatter([col_of[t] for t in dots['taxid1']], [row_of[g] for g in dots['group2']],
               s=dots['source'].map(size), color=DOT_COLOR, linewidths=0, zorder=3)
    ax.set_xlim(-0.6, n_cols - 0.4)
    ax.set_ylim(ys[-1] + 0.6, -0.6)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, linespacing=1.05)
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
    ax.tick_params(axis='x', pad=2 + 2 * strip * figure_style.DENSE_ROW * 72)

    # the margin columns: parasites of each niche, parasite groups
    margin_x = n_cols + 0.6
    headers = [('intracellular', 'intracellular'), ('extracellular', 'extracellular'),
               ('parasite_groups', 'parasite groups')]
    for j, (column, header) in enumerate(headers):
        x = margin_x + 2.4 * j
        # the headers stand like the parasite names, from the top of the strips up
        ax.text(x, -0.6 - 2 * strip - 0.3, header, ha='center', va='bottom',
                rotation=90, color='#555555', clip_on=False)
        for row in ogs.itertuples():
            y = row_of[row.orthologous_group]
            ax.text(x, y, f'{getattr(row, column)}', ha='center', va='center', color='#555555', clip_on=False)

    groups = [Patch(color=color, label=group) for group, color in taxon_colors.items()
              if any(parasites[int(t)]['group'] == group for t in columns)]
    niches = [Patch(color=color, label=niche) for niche, color in NICHE_COLORS.items()]
    sizes = [Line2D([], [], linestyle='', marker='o', color=DOT_COLOR,
                    markersize=size(n) ** 0.5, label=label) for n, label in ((1, '1'), (5, '5'), (20, '20+'))]
    # the keys take the full width: the size key sits beside the niches, as in the GO figure
    legends = figure_style.stack_legends(fig, left=0.01, blocks=[
        ('Taxonomic group', groups), ('Niche', niches)])
    height = fig.get_figheight()
    fig.legend(handles=sizes, title='Parasite proteins', loc='lower left',
               ncol=len(sizes), frameon=False, handlelength=0.9, columnspacing=1.2,
               alignment='left', bbox_to_anchor=(0.705, 0.02 / height),
               title_fontsize=figure_style.FONT)
    # room above the dots for the strips and the parasite names standing on them
    fig.subplots_adjust(left=LEFT, right=0.89, top=1 - 1.2 / height,
                        bottom=(legends + 0.2) / height)
    for extension in ('svg', 'png'):
        fig.savefig(f'{output_stem}.{extension}', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--host', type=int, default=9606, help='host taxid (default: human)')
    parser.add_argument('--min-share', type=float, default=0.5,
                        help='groups reached by a smaller share of the parasites of either niche stay out of the figure')
    parser.add_argument('--table', default='paper/tables/shared_ogs.csv')
    parser.add_argument('--genus-table', default='paper/tables/shared_ogs_genus.csv',
                        help='the groups each genus reaches on its own')
    parser.add_argument('--figure', default='paper/figures/shared_ogs',
                        help='output path without extension')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    predictions = pd.read_parquet(os.path.join(args.data_dir, 'predictions.parquet'))
    ogs, edges = count_ogs(config, predictions, args.host)
    ogs.to_csv(args.table, index=False)
    shared = shared_ogs(ogs, edges, args.min_share)
    print(f"Wrote {args.table}: {len(ogs)} host orthologous groups reached by two or more parasites, "
          f"{len(shared)} by at least {args.min_share:.0%} of each niche")
    genus_ogs = count_genus_ogs(config, edges)
    genus_ogs.to_csv(args.genus_table, index=False)
    counts = genus_ogs.groupby('genus').size()
    print(f"Wrote {args.genus_table}: {len(genus_ogs)} groups reached by one genus alone"
          + (f" ({', '.join(f'{g} {n}' for g, n in counts.items())})" if len(counts) else ''))
    draw(config, shared, edges, args.figure)
    print(f"Wrote {args.figure}.svg and .png")
