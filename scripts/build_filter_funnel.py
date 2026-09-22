'''
Count how many proteins of every species survive each stage of the pipeline, and
draw it. Hosts go proteome -> expressed in an infected tissue -> put by DeepLoc where
a parasite reaches -> in an EggNOG group -> in a group with a high-confidence STRING
link -> in a predicted interaction; parasites go proteome -> secreted/surface ->
group -> linked group -> interaction. Writes paper/tables/funnel.csv (the counts)
and paper/figures/funnel.pdf/.svg (nested bars per species on a log axis, in two
columns).
'''
import argparse
import copy
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils
from pipeline import filters, homology, main
from scripts import figure_style

HOST_STAGES = ['proteome', 'tissue', 'localization', 'eggnog', 'linked', 'interactor']
PARASITE_STAGES = ['proteome', 'localization', 'eggnog', 'linked', 'interactor']
STAGE_LABELS = {
    'proteome': 'STRING proteome',
    'tissue': 'expressed in an infected tissue',
    'localization': 'reachable localization',
    'eggnog': 'in an EggNOG group',
    'linked': 'group has a STRING link ≥ 0.7',
    'interactor': 'in a predicted PPI',
}
# one hue, later stage darker; keyed by stage so both panels agree with the legend
STAGE_COLORS = {
    'proteome': '#e4ecf7',
    'tissue': '#b9d0ec',
    'localization': '#87b0de',
    'eggnog': '#5389c9',
    'linked': '#2a62ab',
    'interactor': '#10315e',
}


def count_stages(config_file, data_dir):
    '''One row per species: the proteins left after every stage.'''
    hosts = utils.read_config(filepath=config_file, field='hosts')
    parasites = utils.read_config(filepath=config_file, field='parasites')
    downloads = os.path.join(data_dir, 'downloads')

    print("Getting proteins...")
    proteins = main.get_proteins(config_file)
    counts = {taxid: {'proteome': len(names)} for taxid, names in proteins.items()}

    print("Applying the tissue filter...")
    tissue_pool = copy.deepcopy(proteins)
    filters.apply_tissue_filter(config_file=config_file, valid_proteins=tissue_pool,
                                cutoff=main.TISSUE_CUTOFF)
    for taxid in hosts:
        counts[taxid]['tissue'] = len(tissue_pool[taxid])

    # eligible_proteins.parquet holds what the secretome, tissue and DeepLoc filters
    # left, so it is the localization stage of both sides
    eligible = pd.read_parquet(os.path.join(data_dir, 'eligible_proteins.parquet'))
    eligible['taxid'] = eligible['taxid'].astype(int)
    for taxid, size in eligible.groupby('taxid').size().items():
        counts[taxid]['localization'] = size

    print("Scanning EggNOG groups...")
    groups = homology.get_eggnog_groups(os.path.join(downloads, '2759_members.tsv.gz'),
                                        eligible['protein'])
    # a protein can sit in several groups
    protein_groups = {}
    for group, members in groups.items():
        for protein in members:
            protein_groups.setdefault(protein, set()).add(group)
    eligible['eggnog'] = eligible['protein'].map(protein_groups)
    for taxid, size in eligible.dropna(subset=['eggnog']).groupby('taxid').size().items():
        counts[taxid]['eggnog'] = size

    print("Scanning COG links...")
    linked = set()
    with utils.read_gzipped_file(os.path.join(downloads, 'COG.links.detailed.v12.0.txt.gz')) as handle:
        next(handle, None)
        for line in handle:
            data = line.decode('utf-8').rstrip().split(' ')
            if data[0] not in groups or data[1] not in groups:
                continue
            cutoff = homology.EVIDENCE_CUTOFF * homology.STRING_SCORE_SCALE
            if int(data[6]) >= cutoff or int(data[7]) >= cutoff:
                linked.update(data[:2])
    has_link = eligible['eggnog'].map(lambda gs: bool(gs & linked), na_action='ignore').fillna(False)
    for taxid, size in eligible[has_link.astype(bool)].groupby('taxid').size().items():
        counts[taxid]['linked'] = size

    predictions = pd.read_parquet(os.path.join(data_dir, 'predictions.parquet'))
    interactors = pd.concat([predictions[['taxid1', 'source']].set_axis(['taxid', 'protein'], axis=1),
                             predictions[['taxid2', 'target']].set_axis(['taxid', 'protein'], axis=1)])
    interactors['taxid'] = interactors['taxid'].astype(int)
    for taxid, size in interactors.drop_duplicates().groupby('taxid').size().items():
        counts[taxid]['interactor'] = size

    group_order = list(utils.read_config(filepath=config_file, field='parasite_groups'))
    rows = []
    for taxid, host in hosts.items():
        rows.append({'taxid': taxid, 'species': host['label'].split(' (')[0], 'side': 'host',
                     'group': 'Host', **counts[taxid]})
    for taxid, parasite in sorted(parasites.items(),
                                  key=lambda item: (group_order.index(item[1]['group']), item[1]['label'])):
        rows.append({'taxid': taxid, 'species': parasite['label'], 'side': 'parasite',
                     'group': parasite['group'], **counts[taxid]})
    table = pd.DataFrame(rows, columns=['taxid', 'species', 'side', 'group'] + HOST_STAGES)
    # a parasite with no prediction has no interactor row; the tissue stage is host-only
    table['interactor'] = table['interactor'].fillna(0)
    for stage in HOST_STAGES:
        table[stage] = table[stage].astype('Int64')

    return table


def draw_side(ax, rows, stages):
    '''Nested bars of the proteins each stage keeps, one row per species, on a log axis.'''
    y = range(len(rows))
    for stage in stages:
        # the log axis starts at 1, so a bar starts there too; a stage with nothing left
        # draws no bar
        ax.barh(y, rows[stage].astype(float).clip(lower=1) - 1, left=1,
                color=STAGE_COLORS[stage], height=0.72, label=STAGE_LABELS[stage])
    ax.set_yticks(list(y))
    ax.set_yticklabels([figure_style.short_name(s) for s in rows['species']], style='italic')
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xscale('log')
    ax.set_xlim(1, 1e5)
    ax.set_xticks([1, 10, 100, 1000, 10000, 100000])
    ax.set_xticklabels(['1', '', '100', '', '10,000', ''])
    ax.minorticks_off()
    ax.tick_params(axis='y', length=0)
    ax.grid(axis='x', color='#e6e6e6', linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_color('#999999')
    for i, (proteome, interactors) in enumerate(zip(rows['proteome'], rows['interactor'])):
        ax.text(1.3e5, i, f'{interactors:,} / {proteome:,}', va='center', color='#555555')


def draw(table, group_colors, output_stem):
    """Two columns: the hosts and the largest parasite group on the left, the rest on the
    right, on a shared grid so every row is the same height."""
    figure_style.apply()
    hosts = table[table['side'] == 'host'].reset_index(drop=True)
    parasites = table[table['side'] == 'parasite']
    largest = parasites['group'].value_counts().idxmax()
    left = parasites[parasites['group'] == largest].reset_index(drop=True)
    right = parasites[parasites['group'] != largest].reset_index(drop=True)

    # one grid row per species; two rows of gap between the panels of the left column
    gap = 2
    n_rows = max(len(hosts) + gap + len(left), len(right))
    fig = plt.figure(figsize=(figure_style.WIDTH, figure_style.ROW * n_rows + 1.9))
    # the gap between the columns holds the left column's counts and the right one's names
    grid = fig.add_gridspec(n_rows, 2, wspace=1.5,
                            top=1 - 0.95 / fig.get_figheight(), bottom=1.05 / fig.get_figheight(),
                            left=0.175, right=0.87)
    top = fig.add_subplot(grid[:len(hosts), 0])
    bottom_left = fig.add_subplot(grid[len(hosts) + gap:len(hosts) + gap + len(left), 0])
    bottom_right = fig.add_subplot(grid[:len(right), 1])

    draw_side(top, hosts, HOST_STAGES)
    draw_side(bottom_left, left, PARASITE_STAGES)
    draw_side(bottom_right, right, PARASITE_STAGES)
    top.set_title('Hosts', loc='left', fontweight='bold')
    bottom_left.set_title('Parasites', loc='left', fontweight='bold')
    bottom_right.set_title('Parasites (continued)', loc='left', fontweight='bold')
    for ax in (bottom_left, bottom_right):
        ax.set_xlabel('Proteins left after each stage')

    # a strip of the taxonomic group colour down the left of its rows, as in the app
    for ax, rows in ((bottom_left, left), (bottom_right, right)):
        for group, members in rows.groupby('group', sort=False):
            first, last = members.index[0], members.index[-1]
            ax.add_patch(Rectangle((-0.82, first - 0.4), 0.025, last - first + 0.8,
                                   color=group_colors[group], linewidth=0,
                                   transform=ax.get_yaxis_transform(), clip_on=False))

    handles, labels = top.get_legend_handles_labels()
    fig.legend(handles, labels, title='Stage', loc='upper center', ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 0.995), handlelength=1.2, columnspacing=1.5,
               alignment='left', title_fontsize=figure_style.FONT)
    figure_style.stack_legends(fig, left=0.175, blocks=[('Taxonomic group', [Patch(color=color, label=group)
                                                          for group, color in group_colors.items()])])
    for extension in ('pdf', 'svg'):
        fig.savefig(f'{output_stem}.{extension}')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--table', default='paper/tables/funnel.csv')
    parser.add_argument('--figure', default='paper/figures/funnel', help='output path without extension')
    args = parser.parse_args()

    table = count_stages(config_file=args.config, data_dir=args.data_dir)
    table.to_csv(args.table, index=False)
    print(f"Wrote {args.table}")
    draw(table, utils.read_config(filepath=args.config, field='parasite_groups'), args.figure)
    print(f"Wrote {args.figure}.pdf and .svg")
