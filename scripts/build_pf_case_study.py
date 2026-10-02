'''
One parasite in one host, narrowed as a funnel: the predicted links at or above a
confidence score, the links of one tissue the parasite infects, and from those two
branches, the links whose host protein is enhanced in one cell type of the tissue and
those whose host protein is found in none of the parasite's other tissues. Draws the
network of each stage (A-D) the way the app draws it, the proteins and links going on to
a later stage in their species colour and the rest muted, and (E) the host processes
most enriched at each stage. The stages are those of scripts/explore_pf_case_study.py. Writes paper/figures/pf_case_study.svg/.png
and the panels' data to paper/tables/pf_case_study_*.csv.
'''
import argparse
import os
import sys
import textwrap

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import networkx as nx
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import utils
from scripts import figure_style
from scripts.build_case_study import HEAT_SCALE, family_name
from scripts.explore_pf_case_study import enhanced, load, tissue_proteins

# the network as the app draws it (app/pages/1_Predicted_Host-Parasite_PPIs.py): a node
# filled with a wash of its species colour inside a border of the colour itself, and
# neutral edges
NODE_FILL_TINT = 0.74
EDGE_COLOR = '#c7cfd9'
MUTED_COLOR = '#c8ced6'
LABEL_FONT_COLOR = '#2f3a46'
# most characters of a parasite protein's description on its node, and per line of it
NODE_LABEL_CHARS = 30
NODE_LABEL_WRAP = 16
# units between the rows of components in the labelled networks, a node row being 0.8
ROW_PITCH = 2.4
# the stages each stage goes on to: the tissue branches into the cell type and the links
# found only in the tissue
CHILDREN = {0: [1], 1: [2, 3]}
# a stage with more proteins than this is drawn unlabelled
MAX_LABELLED = 30
# terms per stage, after dropping those naming the same host proteins as a better one
TERMS_PER_STAGE = 3
FDR = 0.05
# the colour scale of B stops here, so the smaller stages are not all one pale shade
FDR_CAP = 10
# most characters of a term name before it is cut to an ellipsis
TERM_CHARS = 44


def stages(config, data_dir, parasite, host, tissue, cell_type, min_score, min_ntpm=None):
    '''
    The links left at each stage, as [(label, edges)], and the GO annotation. The cell
    type stage keeps the host proteins enhanced in the cell type, or with min_ntpm those
    expressed above it there.
    '''
    edges, annotation, go_df = load(config, data_dir, parasite, host)
    edges = edges.assign(weight=edges['weight'].astype(float))
    edges = edges[edges['weight'] >= min_score]
    in_tissue = tissue_proteins(config, annotation, parasite)
    others = set().union(*(proteins for t, proteins in in_tissue.items() if t != tissue))
    tissue_edges = edges[edges['target'].isin(in_tissue[tissue])]
    if min_ntpm is None:
        in_cell = enhanced(annotation, tissue)[cell_type]
        cell_label = f'Enhanced in\n{cell_type}'
    else:
        cells = annotation[(annotation['Tissue'] == tissue) & (annotation['Cell type'] == cell_type)]
        in_cell = set(cells.loc[cells['nTPM'] > min_ntpm, 'Gene'])
        cell_label = f'{cell_type.capitalize()},\nnTPM > {min_ntpm:g}'
    cell_edges = tissue_edges[tissue_edges['target'].isin(in_cell)]
    specific = tissue_edges[~tissue_edges['target'].isin(others)]

    return [('All\npredicted', edges), (tissue.capitalize(), tissue_edges),
            (cell_label, cell_edges),
            (f'{tissue.capitalize()}-\nonly', specific)], go_df


def enrichment(staged, go_df):
    '''
    The FDR and host protein count of the terms shown, terms by stages. A stage shows its
    best terms, skipping one that names the same host proteins as a term already shown.
    '''
    tests = {label: utils.calculate_enrichment(set(edges['target']), go_df)
             for label, edges in staged}
    shown = []
    for label, test in tests.items():
        significant = test[test['fdr_bh'] < FDR] if not test.empty else test
        seen, picked = set(), 0
        for row in significant.itertuples():
            if picked == TERMS_PER_STAGE:
                break
            if row.nodes in seen:
                continue
            seen.add(row.nodes)
            picked += 1
            if row.go_term not in shown:
                shown.append(row.go_term)
    frames = []
    for label, test in tests.items():
        rows = test[test['go_term'].isin(shown) & (test['fdr_bh'] < FDR)] if not test.empty else test
        frames.append(pd.DataFrame({'term': rows['go_term'], 'stage': label,
                                    'fdr': rows['fdr_bh'], 'proteins': rows['A']}))

    return pd.concat(frames), shown


def strip(ax):
    for side in ax.spines.values():
        side.set_visible(False)


def draw_enrichment(fig, ax, staged, terms, shown):
    '''B: the host processes most enriched at each stage.'''
    rows = {t: i for i, t in enumerate(shown)}
    columns = {label: i for i, (label, _) in enumerate(staged)}
    cmap = LinearSegmentedColormap.from_list('fdr', HEAT_SCALE[1:])
    dots = ax.scatter(terms['stage'].map(columns), terms['term'].map(rows),
                      s=6 * np.sqrt(terms['proteins']) + 4, c=-np.log10(terms['fdr']),
                      cmap=cmap, vmin=-np.log10(FDR), vmax=FDR_CAP, linewidths=0.3,
                      edgecolors='#555555')
    ax.set_yticks(range(len(shown)))
    ax.set_yticklabels([t if len(t) <= TERM_CHARS else t[:TERM_CHARS - 1].rstrip() + '\u2026'
                        for t in shown])
    ax.set_ylim(len(shown) - 0.5, -0.5)
    ax.set_xticks(range(len(staged)))
    ax.set_xticklabels('ABCD'[:len(staged)])
    ax.set_xlabel('stage (network panel)')
    ax.set_xlim(-0.6, len(staged) - 0.4)
    ax.tick_params(length=0)
    ax.grid(color='#eeeeee', linewidth=0.5)
    ax.set_axisbelow(True)
    strip(ax)
    ax.set_title('E  Host processes enriched at each stage', loc='left', pad=4)
    bar = fig.colorbar(dots, ax=ax, fraction=0.04, pad=0.02, anchor=(0, 1), shrink=0.9,
                       extend='max')
    bar.set_label('−log$_{10}$ FDR')
    bar.outline.set_visible(False)
    bar.ax.tick_params(length=2)
    handles = [Line2D([], [], marker='o', linestyle='', markersize=np.sqrt(6 * np.sqrt(n) + 4),
                      markerfacecolor='white', markeredgecolor='#555555', markeredgewidth=0.3,
                      label=str(n)) for n in (2, 10, 100)]
    ax.legend(handles=handles, title='host proteins', loc='upper left',
              bbox_to_anchor=(1.3, 1.05), frameon=False, handletextpad=0.2, labelspacing=0.3)


def tint(color, amount):
    '''A colour mixed toward white, as the app fills its nodes.'''
    channels = [int(color[i:i + 2], 16) for i in (1, 3, 5)]

    return '#%02x%02x%02x' % tuple(round(c + (255 - c) * amount) for c in channels)


def node_labels(edges, descriptions):
    '''
    The name on each node: a host protein's gene name, and a parasite protein's
    description, or its gene name where STRING gives it none.
    '''
    labels = dict(zip(edges['target'], edges['target_name']))
    for protein, name in zip(edges['source'], edges['source_name']):
        description = family_name([descriptions.get(protein, '')])
        description = description.replace(', putative', '')
        labels[protein] = name if description == 'uncharacterized' else (
            textwrap.fill(description if len(description) <= NODE_LABEL_CHARS
                          else description[:NODE_LABEL_CHARS - 1].rstrip() + '\u2026',
                          NODE_LABEL_WRAP, break_on_hyphens=False, break_long_words=False))

    return labels


def packed_layout(graph, aspect, spread, gap, seed=3):
    '''
    Each connected component laid out on its own, sized by its number of nodes, and the
    components packed in rows into a box of the given aspect (width over height), the
    largest first, as the app's physics keeps them in view instead of letting them drift.
    '''
    pieces = []
    for nodes in sorted(nx.connected_components(graph), key=lambda c: (-len(c), sorted(c))):
        sub = graph.subgraph(nodes)
        position = (nx.spring_layout(sub, seed=seed, iterations=200) if len(nodes) > 1
                    else {next(iter(nodes)): np.zeros(2)})
        xy = np.array(list(position.values()), dtype=float)
        xy -= xy.mean(axis=0)
        radius = spread * np.sqrt(len(nodes))
        reach = np.abs(xy).max() or 1
        pieces.append((list(position), xy * radius / reach if len(nodes) > 1 else xy, radius))
    area = sum((2 * r + gap) ** 2 for _, _, r in pieces)
    row_width = np.sqrt(area * aspect)
    placed, x, y, row_height = {}, 0.0, 0.0, 0.0
    for nodes, xy, radius in pieces:
        size = 2 * radius + gap
        if x > 0 and x + size > row_width:
            x, y, row_height = 0.0, y - row_height, 0.0
        for node, (dx, dy) in zip(nodes, xy):
            placed[node] = (x + size / 2 + dx, y - size / 2 + dy)
        x += size
        row_height = max(row_height, size)

    return placed


def two_row_layout(edges, row_width):
    '''
    Each connected component as two rows, its parasite proteins over its host proteins a
    unit apart, ordered so the links cross less; the components packed into rows of at
    most row_width units, each in the first row with room, the largest first. Returns the
    positions and the number of rows.
    '''
    graph = nx.Graph()
    graph.add_edges_from(zip(edges['source'], edges['target']))
    parasites = set(edges['source'])
    placed, filled = {}, []
    for nodes in sorted(nx.connected_components(graph), key=lambda c: (-len(c), sorted(c))):
        top = sorted(n for n in nodes if n in parasites)
        rank = {n: i for i, n in enumerate(top)}
        bottom = sorted((n for n in nodes if n not in parasites),
                        key=lambda h: (np.mean([rank[p] for p in graph[h]]), h))
        size = max(len(top), len(bottom))
        row = next((i for i, used in enumerate(filled) if used + size <= row_width), len(filled))
        if row == len(filled):
            filled.append(0)
        x, y = filled[row], -ROW_PITCH * row
        for nodes_in_row, offset in ((top, 0.0), (bottom, -0.8)):
            for i, node in enumerate(nodes_in_row):
                placed[node] = (x + (i + 0.5) * size / len(nodes_in_row), y + offset)
        filled[row] += size

    return placed, len(filled)


def draw_network(ax, edges, kept, labels, title, labelled, aspect, row_width=None):
    '''
    One stage's network, parasite proteins as diamonds and host proteins as dots, the
    links in kept (those going on to the next stage) and their proteins in colour and the
    rest muted.
    '''
    graph = nx.Graph()
    graph.add_edges_from(zip(edges['source'], edges['target']))
    if labelled:
        position, rows = two_row_layout(edges, row_width)
    else:
        position = packed_layout(graph, aspect, spread=0.3, gap=0.25)
    colors = dict(zip(edges['source'], edges['source_color']))
    colors.update(zip(edges['target'], edges['target_color']))
    kept_pairs = set(zip(kept['source'], kept['target']))
    kept_nodes = set(kept['source']) | set(kept['target'])
    width = 0.3 if not labelled else 0.6
    for row in sorted(edges[['source', 'target', 'weight']].itertuples(index=False),
                      key=lambda r: (r.source, r.target) in kept_pairs):
        on = (row.source, row.target) in kept_pairs
        (x0, y0), (x1, y1) = position[row.source], position[row.target]
        ax.plot([x0, x1], [y0, y1], color=EDGE_COLOR if on else '#e6e9ed',
                linewidth=width * (1 + 2 * (row.weight - 0.7) / 0.3), zorder=1 + on,
                solid_capstyle='round')
    size = 30 if labelled else 9
    for side, marker in (('source', 'D'), ('target', 'o')):
        nodes = sorted(edges[side].unique(), key=lambda n: n in kept_nodes)
        on = [n in kept_nodes for n in nodes]
        ax.scatter([position[n][0] for n in nodes], [position[n][1] for n in nodes],
                   marker=marker, s=size * (0.8 if marker == 'D' else 1),
                   c=[tint(colors[n], NODE_FILL_TINT) if o else '#f1f3f5' for n, o in zip(nodes, on)],
                   edgecolors=[colors[n] if o else MUTED_COLOR for n, o in zip(nodes, on)],
                   linewidths=0.8 if labelled else 0.4, zorder=3)
    if labelled:
        parasites = set(edges['source'])
        for node, (x, y) in position.items():
            above = node in parasites
            ax.annotate(labels[node], (x, y), xytext=(0, 6 if above else -6),
                        textcoords='offset points', ha='center', va='bottom' if above else 'top',
                        color=LABEL_FONT_COLOR, zorder=4, linespacing=1.0)
    if labelled:
        ax.set_xlim(-0.05, row_width + 0.05)
        ax.set_ylim(-ROW_PITCH * (rows - 1) - 1.3, 1.15)
    else:
        ax.set_aspect('equal', adjustable='datalim')
        ax.margins(0.02)
    ax.axis('off')
    ax.set_title(title, loc='left', pad=13, fontsize=figure_style.TITLE)
    ax.text(0, 1.0, f"{len(edges)} links, {edges['source'].nunique()} parasite, "
            f"{edges['target'].nunique()} host proteins", transform=ax.transAxes,
            ha='left', va='bottom', color='#555555')


def draw(staged, terms, shown, labels, parasite_color, host_color, min_score, title,
         output_stem):
    figure_style.apply()
    width = figure_style.WIDTH
    top_h, bottom_h = 2.2, 2.4
    go_h = figure_style.ROW * len(shown) + 0.1
    gaps = [0.6, 0.6, 0.55, 0.65]
    height = top_h + bottom_h + go_h + sum(gaps) + 0.15
    fig = plt.figure(figsize=(width, height))

    def axes(left, top, w, h):
        return fig.add_axes([left / width, 1 - (top + h) / height, w / width, h / height])

    second = gaps[0] + top_h + gaps[1]
    panels = [axes(0.05, gaps[0], 3.4, top_h), axes(3.65, gaps[0], 3.4, top_h),
              axes(0.05, second, 4.85, bottom_h), axes(5.05, second, 2.0, bottom_h)]
    for i, (ax, (label, edges)) in enumerate(zip(panels, staged)):
        later = [staged[c][1] for c in CHILDREN.get(i, [])]
        kept = pd.concat(later) if later else edges
        bounds = ax.get_position()
        draw_network(ax, edges, kept, labels,
                     f'{"ABCD"[i]}  {label.replace(chr(10), " ").replace("- ", "-")}',
                     labelled=edges['source'].nunique() + edges['target'].nunique() <= MAX_LABELLED,
                     aspect=bounds.width * width / (bounds.height * height),
                     row_width=max(1, round(bounds.width * width)))
    handles = [Line2D([], [], marker='D', linestyle='', markersize=4.5,
                      markerfacecolor=tint(parasite_color, NODE_FILL_TINT),
                      markeredgecolor=parasite_color, label='parasite protein'),
               Line2D([], [], marker='o', linestyle='', markersize=5,
                      markerfacecolor=tint(host_color, NODE_FILL_TINT),
                      markeredgecolor=host_color, label='host protein'),
               Line2D([], [], marker='o', linestyle='', markersize=5, markerfacecolor='#f1f3f5',
                      markeredgecolor=MUTED_COLOR, label='left out at the later stages')]
    fig.legend(handles=handles, loc='upper left', ncol=3, frameon=False, handletextpad=0.3,
               bbox_to_anchor=(0.05 / width, 1 - (second + bottom_h + 0.08) / height))

    go = axes(2.65, second + bottom_h + gaps[2], 2.4, go_h)
    draw_enrichment(fig, go, staged, terms, shown)
    fig.suptitle(f'{title}, confidence score \u2265 {min_score:g}', x=0.02, ha='left',
                 y=1 - 0.12 / height, style='italic')
    for extension in ('svg', 'png'):
        fig.savefig(f'{output_stem}.{extension}', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--parasite', type=int, default=36329,
                        help='parasite taxid (default: Plasmodium falciparum)')
    parser.add_argument('--host', type=int, default=9606, help='host taxid (default: human)')
    parser.add_argument('--tissue', default='liver')
    parser.add_argument('--cell-type', default='hepatocytes')
    parser.add_argument('--min-ntpm', type=float,
                        help='keep the host proteins expressed above this nTPM in the cell '
                             'type, instead of those enhanced in it')
    parser.add_argument('--min-score', type=float, default=0.7,
                        help='links below this confidence score are left out')
    parser.add_argument('--tables', default='paper/tables/pf_case_study',
                        help='stem of the tables written')
    parser.add_argument('--figure', default='paper/figures/pf_case_study',
                        help='output path without extension')
    args = parser.parse_args()

    config = utils.read_config(filepath=args.config)
    staged, go_df = stages(config, args.data_dir, args.parasite, args.host, args.tissue,
                           args.cell_type, args.min_score, args.min_ntpm)
    terms, shown = enrichment(staged, go_df)
    descriptions = pd.read_parquet(os.path.join(args.data_dir, 'protein_annotations.parquet')
                                   ).set_index('protein')['description']
    labels = node_labels(staged[0][1], descriptions)
    title = (f"{config['parasites'][args.parasite]['label']} in "
             f"{config['hosts'][args.host]['label'].split(' (')[0]}")
    draw(staged, terms, shown, labels, config['parasites'][args.parasite]['color'],
         staged[0][1]['target_color'].iloc[0], args.min_score, title, args.figure)

    names = [label.replace('\n', ' ').replace('- ', '-') for label, _ in staged]
    counts = pd.DataFrame([{'stage': n, 'links': len(e), 'parasite proteins': e['source'].nunique(),
                            'host proteins': e['target'].nunique()} for n, (_, e) in zip(names, staged)])
    counts.to_csv(f'{args.tables}_funnel.csv', index=False)
    terms.assign(stage=terms['stage'].map(dict(zip([l for l, _ in staged], names)))
                 ).to_csv(f'{args.tables}_go.csv', index=False)
    pd.concat([edges.assign(stage=name) for name, (_, edges) in zip(names[2:], staged[2:])]
              ).assign(parasite_label=lambda d: d['source'].map(labels))[
        ['stage', 'source', 'source_name', 'parasite_label', 'target', 'target_name', 'weight']
    ].to_csv(f'{args.tables}_links.csv', index=False)
    print(counts.to_string(index=False))
    print(f"Wrote {args.figure}.svg and .png and {args.tables}_*.csv")
