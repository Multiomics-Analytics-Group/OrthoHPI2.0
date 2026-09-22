'''
Compare the hosts of every parasite predicted against more than one. An interaction is
compared between hosts as the pair of orthologous groups it was transferred from, so a
link counts as shared when both hosts carry a protein of the host group interacting with
the parasite's. A link missing from a host is explained with the first filter its host
family fails there: the family is absent from the host proteome, no member is expressed
in a tissue the parasite infects, no expressed member is in a location the parasite
reaches, or it is available and STRING transferred nothing. Writes
paper/tables/multi_host.csv and paper/figures/multi_host.pdf/.svg.
'''
import argparse
import copy
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils
from pipeline import filters, main
from scripts import figure_style

REASONS = ['absent', 'not expressed', 'out of reach', 'not transferred']
REASON_LABELS = {
    'absent': 'family absent from the host',
    'not expressed': 'not expressed in an infected tissue',
    'out of reach': 'not in a reachable location',
    'not transferred': 'available, not transferred',
}
REASON_COLORS = {'absent': '#3b3b3b', 'not expressed': '#7a7a7a',
                 'out of reach': '#b4b4b4', 'not transferred': '#e0e0e0'}
ALL_HOSTS_COLOR = '#3b3b3b'
SOME_HOSTS_COLOR = '#b4b4b4'


def host_pools(config_file, data_dir):
    '''Per host protein: the infected-tissue labels it is expressed in, and the niches
    that reach it, straight from the TISSUES and DeepLoc files the pipeline reads.'''
    proteins = main.get_proteins(config_file)
    print("Applying the tissue filter...")
    tissues = filters.apply_tissue_filter(config_file=config_file,
                                          valid_proteins=copy.deepcopy(proteins),
                                          cutoff=main.TISSUE_CUTOFF)
    print("Applying the DeepLoc filter...")
    reachable = filters.apply_deeploc_filter(
        config_file=config_file, valid_proteins=copy.deepcopy(proteins),
        deeploc_dir=os.path.join(data_dir, main.DEEPLOC_ACCURATE_DIR),
        niche_cutoffs=main.DEEPLOC_NICHE_CUTOFFS, default_niche=main.DEEPLOC_DEFAULT_NICHE)

    return tissues, reachable


def family_members(data_dir):
    '''{(host group, host taxid): its proteins in that host}.'''
    orthologs = pd.read_parquet(os.path.join(data_dir, 'host_orthologs.parquet'))
    members = {}
    for row in orthologs.itertuples():
        members.setdefault((row.group, str(row.taxid)), set()).update(row.proteins.split(','))

    return members


def group_links(edges):
    '''The predictions of one parasite as unordered group pairs, with the hosts carrying
    each and the host group of the pair.'''
    edges = edges.copy()
    edges['pair'] = [tuple(sorted(pair)) for pair in zip(edges['group1'], edges['group2'])]

    return edges.groupby('pair').agg(hosts=('taxid2', frozenset),
                                     host_group=('group2', lambda g: g.mode().iloc[0]))


def why_missing(link, host, parasite, members, tissues, reachable):
    family = members.get((link.host_group, host), set())
    if not family:
        return 'absent'
    expressed = {p for p in family if parasite['tissue_labels'] & set(tissues.get(p, []))}
    if not expressed:
        return 'not expressed'
    if not expressed & reachable[parasite['niche']]:
        return 'out of reach'
    return 'not transferred'


def compare(config_file, data_dir):
    '''One row per (parasite, host) of the multi-host parasites.'''
    config = utils.read_config(filepath=config_file)
    hosts = {str(t): h['label'].split('(')[1].rstrip(')') for t, h in config['hosts'].items()}
    predictions = pd.read_parquet(os.path.join(data_dir, 'predictions.parquet'))
    tissues, reachable = host_pools(config_file, data_dir)
    members = family_members(data_dir)

    rows = []
    for taxid, parasite in config['parasites'].items():
        if len(parasite['hosts']) < 2:
            continue
        parasite = dict(parasite, tissue_labels={config['tissues'][t] for t in parasite['tissues']})
        parasite_hosts = [str(h) for h in parasite['hosts']]
        edges = predictions[predictions['taxid1'] == str(taxid)]
        links = group_links(edges)
        n_hosts = links['hosts'].map(len)
        for host in parasite_hosts:
            in_host = links['hosts'].map(lambda h: host in h)
            reasons = [why_missing(link, host, parasite, members, tissues, reachable)
                       for link in links[~in_host].itertuples()]
            rows.append({
                'parasite': parasite['label'], 'host': hosts[host],
                'predicted PPIs': int((edges['taxid2'] == host).sum()),
                'links': int(in_host.sum()),
                'links in all hosts': int((n_hosts == len(parasite_hosts)).sum()),
                'links in some hosts': int((in_host & (n_hosts > 1) & (n_hosts < len(parasite_hosts))).sum()),
                'links only here': int((in_host & (n_hosts == 1)).sum()),
                **{f'missing: {reason}': reasons.count(reason) for reason in REASONS}})
        print(f"{parasite['label']}: {len(links)} links, "
              f"{int((n_hosts == len(parasite_hosts)).sum())} in every host")

    return pd.DataFrame(rows)


def draw(table, config, output_stem):
    figure_style.apply()
    host_colors = {h['label'].split('(')[1].rstrip(')'): h['color'] for h in config['hosts'].values()}
    parasites = list(dict.fromkeys(table['parasite']))
    fig, (left, right) = plt.subplots(1, 2, figsize=(figure_style.WIDTH, 0.5 * len(parasites) + 1.9),
                                      gridspec_kw={'width_ratios': [1, 1], 'wspace': 0.22})

    # A: the group-pair links of each parasite, by the hosts carrying them
    for i, parasite in enumerate(parasites):
        rows = table[table['parasite'] == parasite]
        first = rows.iloc[0]
        x = 0
        for value, color in ((first['links in all hosts'], ALL_HOSTS_COLOR),
                             (rows['links in some hosts'].sum() // 2 if len(rows) > 2 else 0, SOME_HOSTS_COLOR)):
            # every link in some hosts of a three-host parasite is counted by two hosts
            left.barh(i, value, left=x, color=color, height=0.7)
            x += value
        for row in rows.to_dict('records'):
            value = row['links only here']
            left.barh(i, value, left=x, color=host_colors[row['host']], height=0.7)
            x += value
        left.text(x + 3, i, f'{x:,}', va='center', color='#555555')

    # B: the links each host lacks while another host of the parasite has them, by cause,
    # one thin bar per host on the parasite's row
    ticks, labels = [], []
    for i, parasite in enumerate(parasites):
        rows = table[table['parasite'] == parasite].to_dict('records')
        height = 0.7 / len(rows)
        for j, row in enumerate(rows):
            y = i - 0.35 + height * (j + 0.5)
            x = 0
            for reason in REASONS:
                value = row[f'missing: {reason}']
                right.barh(y, value, left=x, color=REASON_COLORS[reason], height=height * 0.9)
                x += value
            right.text(x + 3, y, f'{x:,}', va='center', color='#555555')
            ticks.append(y)
            labels.append(row['host'])

    left.set_yticks(range(len(parasites)))
    left.set_yticklabels([figure_style.short_name(p) for p in parasites], style='italic')
    right.set_yticks(ticks)
    right.set_yticklabels(labels)
    for ax in (left, right):
        ax.set_ylim(len(parasites) - 0.4, -0.6)
    left.set_title('A  Links of each parasite, by host', loc='left', fontweight='bold')
    right.set_title('B  Missing from a host, by cause', loc='left', fontweight='bold')
    for ax in (left, right):
        ax.tick_params(axis='y', length=0)
        ax.grid(axis='x', color='#e6e6e6', linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ('top', 'right', 'left'):
            ax.spines[side].set_visible(False)
        ax.spines['bottom'].set_color('#999999')
        ax.set_xlabel('Group-pair links')
        ax.margins(x=0.2)

    host_handles = [Patch(color=ALL_HOSTS_COLOR, label='in every host'),
                    Patch(color=SOME_HOSTS_COLOR, label='in two of three hosts')]
    host_handles += [Patch(color=color, label=f'only in {host}') for host, color in host_colors.items()
                     if host in set(table['host'])]
    left.legend(handles=host_handles, loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=2,
                frameon=False, handlelength=0.9)
    right.legend(handles=[Patch(color=REASON_COLORS[r], label=REASON_LABELS[r]) for r in REASONS],
                 loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=1,
                 frameon=False, handlelength=0.9)
    fig.subplots_adjust(left=0.15, right=0.95, top=1 - 0.35 / fig.get_figheight(),
                        bottom=1.35 / fig.get_figheight())
    for extension in ('pdf', 'svg'):
        fig.savefig(f'{output_stem}.{extension}')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--table', default='paper/tables/multi_host.csv')
    parser.add_argument('--figure', default='paper/figures/multi_host', help='output path without extension')
    args = parser.parse_args()

    table = compare(config_file=args.config, data_dir=args.data_dir)
    table.to_csv(args.table, index=False)
    print(f"Wrote {args.table}")
    draw(table, utils.read_config(filepath=args.config), args.figure)
    print(f"Wrote {args.figure}.pdf and .svg")
