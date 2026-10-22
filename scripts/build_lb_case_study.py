'''
One parasite in its two hosts. An interaction is compared between the hosts as the pair
of orthologous groups it was transferred from, as in scripts/build_multi_host_figure.py,
and a link one host lacks is explained by the first filter its host group fails there.
Draws (A) the app's body figure of each host, its infected organs shaded by the
interactions reaching them, (B) the links of each host, shared or found only there, the
latter by why the other host lacks them, and (C, D) for chosen parasite orthologous groups,
the host groups each host links them to, coloured by confidence score, with the cause
marked where a host lacks the link. The body figures are rendered with Inkscape. Writes paper/figures/lb_case_study.svg/.png and the panels' data to
paper/tables/lb_case_study_*.csv.
'''
import argparse
import os
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'app'))

import body_figure
import utils
from scripts import figure_style
from scripts.build_case_study import HEAT_SCALE, og_name
from scripts.build_multi_host_figure import (REASON_COLORS, REASON_EDGE, REASON_LABELS,
                                             group_members, group_links, host_pools,
                                             why_missing)

SHARED_COLOR = '#3b3b3b'
# the mark in a matrix cell a host lacks, by cause
REASON_MARKS = {'absent': '–', 'not expressed': 'T', 'out of reach': 'L',
                'not transferred': '?'}
REASON_SHORT = {'absent': 'group absent', 'not expressed': 'TISSUES filter',
                'out of reach': 'DeepLoc filter', 'not transferred': 'not transferred'}
INKSCAPE = '/Applications/Inkscape.app/Contents/MacOS/inkscape'
# host gene names listed on a column before the rest are counted
LABEL_NAMES = 2
# inches a character of an 8 pt column label takes
CHAR_WIDTH = 0.068


def host_name(config, taxid):
    return config['hosts'][int(taxid)]['label'].split('(')[1].rstrip(')')


def compare(config_file, config, data_dir, parasite_taxid):
    '''
    The parasite's predictions, and its group-pair links with the hosts carrying each
    and, for the host lacking it, why.
    '''
    parasite = config['parasites'][parasite_taxid]
    hosts = [str(h) for h in parasite['hosts']]
    if len(hosts) != 2:
        sys.exit(f"{parasite['label']} has {len(hosts)} hosts; this figure compares two")
    parasite = dict(parasite, tissue_labels={config['tissues'][t] for t in parasite['tissues']})
    predictions = pd.read_parquet(os.path.join(data_dir, 'predictions.parquet'))
    edges = predictions[predictions['taxid1'] == str(parasite_taxid)].copy()
    edges['weight'] = edges['weight'].astype(float)
    links = group_links(edges)
    tissues, reachable = host_pools(config_file, data_dir)
    members = group_members(data_dir)
    links['missing'] = [None if len(link.hosts) == 2 else
                        why_missing(link, next(h for h in hosts if h not in link.hosts),
                                    parasite, members, tissues, reachable)
                        for link in links.itertuples()]

    return edges, links, hosts


def short_names(names):
    '''Gene names differing only in a final letter, or a final digit after a letter, as
    one: ATG4A, ATG4B and ATG4D as ATG4A/B/D.'''
    groups = {}
    for name in sorted(set(names)):
        stem = name[:-1]
        key = stem if len(name) > 2 and stem[-1:].isalpha() or name[-1:].isalpha() else name
        groups.setdefault(key, []).append(name)
    parts = [group[0] if len(group) == 1 else group[0] + '/' + '/'.join(n[-1] for n in group[1:])
             for group in groups.values()]
    if len(parts) > LABEL_NAMES:
        parts = parts[:LABEL_NAMES] + [f'+{len(parts) - LABEL_NAMES}']

    return ', '.join(parts)


def og_matrix(edges, links, hosts, group, max_columns):
    '''
    The host orthologous groups a parasite group links to, hosts by host groups: the
    best confidence score of the link in each host, and the cause where a host lacks it.
    The max_columns groups with the best score are kept.
    '''
    linked = edges[(edges['group1'] == group) | (edges['group2'] == group)]
    host_group = np.where(linked['group1'] == group, linked['group2'], linked['group1'])
    linked = linked.assign(host_group=host_group)
    scores = linked.pivot_table(index='taxid2', columns='host_group', values='weight',
                                aggfunc='max').reindex(hosts)
    order = scores.max().sort_values(ascending=False, kind='stable').index[:max_columns]
    scores = scores[order]
    names = linked.groupby('host_group')['target_name'].agg(
        lambda n: short_names(n.str.upper()))[order]
    pair = {h: tuple(sorted((group, h))) for h in order}
    causes = pd.DataFrame({h: [None if not np.isnan(scores.at[t, h]) else
                               links.at[pair[h], 'missing'] for t in hosts] for h in order},
                          index=hosts)

    return scores, causes, names


def body_images(config, data_dir, edges, hosts, inkscape):
    """
    The app's body figure of each host, its infected organs shaded on one scale by the
    interactions whose host protein is annotated to them, as images; the organ counts of
    each host; and the colour bins.
    """
    figure_tissues = pd.read_parquet(os.path.join(data_dir, 'figure_tissues.parquet'))
    infected = body_figure.infected_organs(config, edges['taxid1'].iloc[0])
    drawn = []
    for host in hosts:
        species = body_figure.get_species(config, int(host))
        svg, organs = body_figure.load_figure(species)
        counts = body_figure.count_interactions(edges[edges['taxid2'] == host], figure_tissues)
        drawn.append((species, svg, {o: n for o, n in counts.items()
                                     if o in organs and o in infected}))
    bounds, _ = body_figure.color_scale({(i, o): n for i, (_, _, counts) in enumerate(drawn)
                                         for o, n in counts.items()})
    images = []
    with tempfile.TemporaryDirectory() as tmp:
        for species, svg, counts in drawn:
            root = ET.fromstring(body_figure.shade_figure(svg, counts, bounds))
            if species == 'human':
                _, y, _, h = (float(v) for v in root.get('viewBox').replace(',', ' ').split())
                x, w = body_figure.HUMAN_FRONT_VIEW
                root.set('viewBox', f'{x} {y} {w} {h}')
            root.attrib.pop('style', None)
            source, target = os.path.join(tmp, f'{species}.svg'), os.path.join(tmp, f'{species}.png')
            with open(source, 'w') as handle:
                handle.write(ET.tostring(root, encoding='unicode'))
            subprocess.run([inkscape, source, '--export-type=png', '--export-dpi=300',
                            '--export-background=white', f'--export-filename={target}'],
                           check=True, capture_output=True)
            images.append(plt.imread(target))

    return images, [counts for _, _, counts in drawn]


def draw_bodies(axes_list, images, counts, hosts, names):
    """A: the body figure of each host, with its organ counts under it."""
    for ax, image, host_counts, host in zip(axes_list, images, counts, hosts):
        ax.imshow(image, interpolation='lanczos')
        ax.axis('off')
        ax.set_title(names[host], fontsize=figure_style.FONT, pad=2)
        ax.text(0.5, -0.04, '\n'.join(f'{o} {n}' for o, n in
                                      sorted(host_counts.items(), key=lambda c: -c[1])),
                transform=ax.transAxes, ha='center', va='top', color='#555555')

def draw_counts(ax, links, hosts, names):
    '''A: each host's links, shared or found only there, by why the other host lacks them.'''
    shared = int((links['hosts'].map(len) == 2).sum())
    for i, host in enumerate(hosts):
        only = links[links['hosts'].map(lambda h: h == frozenset([host]))]
        ax.barh(i, shared, color=SHARED_COLOR, height=0.6)
        x = shared
        for reason in REASON_COLORS:
            value = int((only['missing'] == reason).sum())
            if value:
                ax.barh(i, value, left=x, color=REASON_COLORS[reason], height=0.6,
                        edgecolor=REASON_EDGE, linewidth=0.4)
                x += value
        ax.text(x + 3, i, f'{shared} + {len(only)}',
                va='center', color='#555555')
    ax.set_yticks(range(len(hosts)))
    ax.set_yticklabels([names[h] for h in hosts])
    ax.set_ylim(len(hosts) - 0.4, -0.6)
    ax.set_xlabel('Group-pair links')
    ax.margins(x=0.12)
    ax.tick_params(axis='y', length=0)
    ax.grid(axis='x', color='#e6e6e6', linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_color('#999999')
    ax.set_title('B  Links of each host', loc='left')
    drawn = [r for r in REASON_COLORS if links['missing'].eq(r).any()]
    handles = [Patch(color=SHARED_COLOR, label='in both hosts')]
    handles += [Patch(facecolor=REASON_COLORS[r], edgecolor=REASON_EDGE, linewidth=0.4,
                      label=f'only in this host; in the other, {REASON_SHORT[r]} '
                            f'({REASON_MARKS[r]})')
                for r in drawn]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(-0.2, -0.75), frameon=False,
              handlelength=0.9, labelspacing=0.3)


def draw_matrix(ax, scores, causes, names, hosts, host_names, title, cmap, vmin):
    '''A parasite orthologous group's host groups: hosts as rows, confidence score as colour.'''
    values = scores.to_numpy(dtype=float)
    image = ax.imshow(np.ma.masked_invalid(values), cmap=cmap, vmin=vmin, vmax=1,
                      aspect='auto')
    for i, host in enumerate(hosts):
        for j, column in enumerate(scores.columns):
            cause = causes.at[host, column]
            if cause:
                ax.text(j, i, REASON_MARKS[cause], ha='center', va='center', color='#555555')
    ax.set_xticks(range(len(scores.columns)))
    ax.set_xticklabels(names, rotation=90)
    ax.set_yticks(range(len(hosts)))
    ax.set_yticklabels([host_names[h] for h in hosts])
    ax.set_xticks(np.arange(-0.5, len(scores.columns)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(hosts)), minor=True)
    ax.grid(which='minor', color='#dddddd', linewidth=0.5)
    ax.tick_params(which='both', length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.set_title(title, loc='left')

    return image


def draw(config, edges, links, hosts, matrices, bodies, title, output_stem):
    figure_style.apply()
    width = figure_style.WIDTH
    names = {h: host_name(config, h) for h in hosts}
    cell = figure_style.ROW
    top_h = 1.5
    matrix_h = len(hosts) * cell
    gaps = [0.75, 0.95]
    label_h = CHAR_WIDTH * max(len(n) for m in matrices for n in m[3]) + 0.1
    height = gaps[0] + top_h + gaps[1] + matrix_h + label_h + 0.05
    fig = plt.figure(figsize=(width, height))

    def axes(left, top, w, h):
        return fig.add_axes([left / width, 1 - (top + h) / height, w / width, h / height])

    images, counts = bodies
    # a portrait figure (the human) takes the row's height, a landscape one two thirds of it
    heights = [top_h if i.shape[0] > i.shape[1] else top_h * 2 / 3 for i in images]
    x, body_axes = 0.2, []
    for image, h in zip(images, heights):
        w = h * image.shape[1] / image.shape[0]
        body_axes.append(axes(x, gaps[0] + top_h - h, w, h))
        x += w + 0.25
    draw_bodies(body_axes, images, counts, hosts, names)
    fig.text(0.05 / width, 1 - (gaps[0] - 0.15) / height,
             'A  Predicted interactions per infected organ', fontsize=figure_style.TITLE,
             va='bottom')

    left = 0.55
    draw_counts(axes(3.95, gaps[0] + 0.05, 2.0, 0.5), links, hosts, names)
    cmap = LinearSegmentedColormap.from_list('score', HEAT_SCALE[1:])
    vmin = float(edges['weight'].min())
    top = gaps[0] + top_h + gaps[1]
    x = left
    for i, (letter, (label, scores, causes, column_names)) in enumerate(zip('CD', matrices)):
        ax = axes(x, top, cell * scores.shape[1], matrix_h)
        image = draw_matrix(ax, scores, causes, column_names, hosts, names,
                            f'{letter}  {label}', cmap, vmin)
        if i:
            ax.set_yticklabels([])
        x += cell * scores.shape[1] + 0.3
    bar_ax = axes(3.6, gaps[0] + top_h + 0.03, 1.4, 0.08)
    bar = fig.colorbar(image, cax=bar_ax, orientation='horizontal')
    bar.set_label('confidence score (C, D)')
    bar.outline.set_visible(False)
    bar.ax.tick_params(length=2)
    fig.suptitle(title, x=0.02, ha='left', y=1 - 0.08 / height, style='italic')
    for extension in ('svg', 'png'):
        fig.savefig(f'{output_stem}.{extension}', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--parasite', type=int, default=5660,
                        help='parasite taxid (default: Leishmania braziliensis)')
    parser.add_argument('--groups', nargs='+', default=['KOG1654', 'KOG2556'],
                        help='parasite orthologous groups drawn in C, D (default: ATG8, GP63)')
    parser.add_argument('--max-columns', type=int, default=22,
                        help='most host groups drawn per parasite group, the best-scoring')
    parser.add_argument('--inkscape', default=INKSCAPE, help='renders the body figures')
    parser.add_argument('--tables', default='paper/tables/lb_case_study',
                        help='stem of the tables written')
    parser.add_argument('--figure', default='paper/figures/lb_case_study',
                        help='output path without extension')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    edges, links, hosts = compare(args.config, config, args.data_dir, args.parasite)
    descriptions = pd.read_parquet(os.path.join(args.data_dir, 'protein_annotations.parquet')
                                   ).set_index('protein')['description']
    matrices, rows = [], []
    for group in args.groups:
        proteins = edges.loc[edges['group1'] == group, 'source'].unique()
        if not len(proteins):
            sys.exit(f'{group} has no predictions for parasite {args.parasite}')
        name = og_name([descriptions.get(p, '') for p in proteins])
        scores, causes, column_names = og_matrix(edges, links, hosts, group, args.max_columns)
        n = len(proteins)
        matrices.append((f'{name}\n{n} parasite protein{"s" * (n > 1)}, {group}',
                         scores, causes, column_names))
        for host in hosts:
            for column in scores.columns:
                rows.append({'parasite group': group, 'parasite group name': name,
                             'host group': column, 'host genes': column_names[column],
                             'host': host_name(config, host), 'score': scores.at[host, column],
                             'missing': causes.at[host, column]})
    title = (f"{config['parasites'][args.parasite]['label']} in "
             f"{' and '.join(host_name(config, h) for h in hosts)}")
    bodies = body_images(config, args.data_dir, edges, hosts, args.inkscape)
    draw(config, edges, links, hosts, matrices, bodies, title, args.figure)

    pd.DataFrame(rows).to_csv(f'{args.tables}_ogs.csv', index=False)
    only = links[links['hosts'].map(len) == 1].reset_index()
    host_genes = edges.groupby(['group1', 'group2', 'taxid2'])['target_name'].agg(
        lambda n: ','.join(sorted(set(n))))
    parasite_genes = edges.groupby(['group1', 'group2'])['source_name'].agg(
        lambda n: ','.join(sorted(set(n))))
    specific = pd.DataFrame([{
        'only in': host_name(config, next(iter(row.hosts))),
        'missing from the other host': row.missing,
        'host group': row.host_group,
        'parasite proteins': next(parasite_genes.get(k) for k in (row.pair, row.pair[::-1])
                                  if k in parasite_genes.index),
        'host proteins': next(host_genes.get(k + (next(iter(row.hosts)),))
                              for k in (row.pair, row.pair[::-1])
                              if k + (next(iter(row.hosts)),) in host_genes.index),
        'score': edges.loc[edges['group1'].isin(row.pair) & edges['group2'].isin(row.pair),
                           'weight'].max(),
    } for row in only.itertuples()]).sort_values(['only in', 'missing from the other host', 'score'],
                                                  ascending=[True, True, False])
    specific.to_csv(f'{args.tables}_specific.csv', index=False)
    print(specific.groupby(['only in', 'missing from the other host']).size().to_string())
    print(f"Wrote {args.figure}.svg and .png and {args.tables}_*.csv")
