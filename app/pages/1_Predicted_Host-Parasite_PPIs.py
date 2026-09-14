import sys, os
import json
import re
import textwrap
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import utils
import web_utils
import streamlit as st
from st_aggrid import GridOptionsBuilder, AgGrid
import pandas as pd
import numpy as np
import networkx as nx
from css import style
from pyvis.network import Network
import plotly.express as px
import plotly.graph_objects as pgo
import structure_visualizer as strv
import body_figure
from ppi_network import ppi_network

style.load_css()
web_utils.show_header('Host-parasite network')

df_select = None
networks = []
# named before the column that draws them, so the sections below can say what they show
selected_tissues = []
selected_cell_types = []
selected_surface = []
surface_options = []
selected_rows = []
selected_terms = []
enrichment_table = None
enrichment_view = None
enrichment = None
path = 'data/tmp'
# where the network is written for the download buttons; not in the repository
os.makedirs(path, exist_ok=True)
# characters per line of a node label
LABEL_WRAP_WIDTH = 22
LABEL_FONT_SIZE = 22  # vis.js defaults to 14, too small to read the protein names
LABEL_FONT_COLOR = '#2f3a46'
LABEL_FONT_FACE = "'Source Sans Pro', -apple-system, 'Segoe UI', sans-serif"
# halo behind the labels, so they are not read through the edges
NETWORK_BACKGROUND = '#fbfcfd'
# how far toward white a node's species colour is washed out to fill it, so the label stays
# readable
NODE_FILL_TINT = 0.74
NODE_HIGHLIGHT_TINT = 0.5
# edges are neutral until hovered
EDGE_COLOR = '#c7cfd9'
EDGE_ACCENT_COLOR = '#2b8cbe'
# the enrichment view, where selected processes are picked out of a dimmed network
HIGHLIGHT_COLOR = '#e7298a'
MUTED_COLOR = '#c8ced6'
# significance is a magnitude, so a single hue running light to dark
GO_SEQUENTIAL = ['#bcdcec', '#7fc0dd', '#3f9fca', '#2b8cbe', '#12587d', '#08324a']
GO_AXIS_COLOR = '#8d97a3'
# fewer proteins than this and the cell-type matrix says nothing the table does not
MATRIX_MIN_PROTEINS = 3
MATRIX_MARK_SIZE = 11
# the bars do not grow with the host proteins
BARS_HEIGHT = 420
# joins tissue and cell type, a cell type name being unique only inside its tissue
MATRIX_SEPARATOR = ' | '
GO_GRID_COLOR = '#e6eaef'
# processes the ranked plot shows; beyond this the names stop being readable
GO_TOP_N = 20
# characters per line of a GO term name on the ranked plot
GO_LABEL_WRAP_WIDTH = 42
# an FDR that underflows to 0 would be an infinite -log10
GO_MIN_FDR = 1e-300
# the treemap root, in the data so plotly does not paint it grey
GO_TREEMAP_ROOT_ID = '__all_enriched_processes__'
GO_TREEMAP_ROOT_LABEL = 'All enriched processes'
# the layout spaces the proteins widely; a shorter canvas only shrinks the whole network
NETWORK_HEIGHT = 1000
# height of each structure viewer in the interaction dialog
VIEWER_HEIGHT = 400

# columns of the interactions table; both STRING ids are carried since the host one cannot
# be assembled from the rest
TABLE_COLUMNS = ['source_name', 'source_full_name', 'source', 'source_uniprot',
                 'target_name', 'target_full_name', 'target', 'target_uniprot',
                 'Tissues', 'experimental_evidence_score',
                 'databases_evidence_score', 'weight', 'group1', 'group2']

# `source` is the parasite and `target` the host throughout the pipeline; the names here are
# the ones the rest of the app uses
TABLE_COLUMN_NAMES = {
    'source_name': 'Parasite protein',
    'source_full_name': 'Parasite protein description',
    'source': 'Parasite STRING id',
    'source_uniprot': 'Parasite UniProt id',
    'source_surface': 'Parasite DeepLoc class',
    'taxid2_label': 'Host species',
    'target_name': 'Host protein',
    'target_full_name': 'Host protein description',
    'target': 'Host STRING id',
    'target_uniprot': 'Host UniProt id',
    'target_surface': 'Host DeepLoc class',
    'experimental_evidence_score': 'Experimental evidence',
    'databases_evidence_score': 'Database evidence',
    'weight': 'Confidence score',
    'group1': 'Parasite orthology group',
    'group2': 'Host orthology group'}

# GO_TERM_COLUMN is read back out of the grid when rows are picked
ENRICHMENT_COLUMN_NAMES = {
    'go_term': 'Biological process',
    'n_proteins': 'Proteins of the network',
    'odds_ratio': 'Odds ratio',
    'p_value': 'P value',
    'fdr_bh': 'FDR (BH)'}
GO_TERM_COLUMN = ENRICHMENT_COLUMN_NAMES['go_term']

# filters written as badges above each section, one colour per kind; the score is always one
# of them
FILTER_BADGES = {'tissue': 'blue', 'cell type': 'green', 'localisation': 'violet'}
# names of one kind in a download filename before they are counted instead
FILENAME_MAX_NAMES = 2

# rows to a page of either table
TABLE_PAGE_SIZE = 10

# the two sides of the network, tested separately for enrichment, also when both are shown
HOST, PARASITE = 'Host proteins', 'Parasite proteins'
BOTH = 'Both sides'
SIDE_COLUMN = 'Side'

# named so a click on the body figure can set it on the following run
TISSUE_FILTER_KEY = 'net_tissues'

# named so their state can be read where the predictions are filtered, before they are drawn
# again
SURFACE_FILTER_KEYS = {web_utils.CELL_MEMBRANE: 'net_surface_membrane',
                       web_utils.EXTRACELLULAR: 'net_surface_extracellular',
                       web_utils.CYTOPLASM: 'net_surface_cytoplasm',
                       web_utils.NUCLEUS: 'net_surface_nucleus'}

config = utils.read_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()


# loaded by web_utils, so the pages share one cached copy
load_predictions = web_utils.load_predictions


@st.cache_data(show_spinner=False, max_entries=5)
def get_parasite_tissues(data_dir, parasite, host_taxids):
    '''
    Annotates the predictions of one parasite against the selected host with the tissues
    and cell types its host targets are expressed in.
    '''
    predictions = load_predictions(data_dir)
    tissues = utils.read_parquet_file(input_file=f'{data_dir}/tissues_cell_types.parquet')
    predictions = predictions[(predictions['taxid1_label'] == parasite)
                              & (predictions['taxid2'].isin(host_taxids))]

    return pd.merge(predictions, tissues.rename({'Gene': 'target'}, axis=1), on='target', how='left')


@st.cache_data(show_spinner=False)
def get_parasite_list(data_dir, host_taxids):
    '''
    Parasites with predictions against the selected host, so the parasite selector never
    offers one that yields an empty network.
    '''
    predictions = load_predictions(data_dir)
    predictions = predictions[predictions['taxid2'].isin(host_taxids)]

    return predictions['taxid1_label'].sort_values().unique().tolist()


@st.cache_data(show_spinner=False)
def load_ontology(data_dir):
    return utils.read_parquet_file(input_file=f'{data_dir}/go_ontology.parquet')


def tint(color, amount):
    '''
    Mixes a colour toward white, so the colour a species is identified by can also be
    used as a fill light enough to write on.
    '''
    color = str(color)
    if not color.startswith('#') or len(color) != 7:
        return color
    channels = [int(color[i:i + 2], 16) for i in (1, 3, 5)]

    return '#%02x%02x%02x' % tuple(round(c + (255 - c) * amount) for c in channels)


def node_color(color):
    '''
    The fill, border and hover colours of a node from the one colour its species is
    drawn in: a wash of the colour inside a border of the colour itself, which reads as
    a group at a glance and keeps the protein name on top of it legible.
    '''
    return {'background': tint(color, NODE_FILL_TINT),
            'border': color,
            'highlight': {'background': tint(color, NODE_HIGHLIGHT_TINT), 'border': color},
            'hover': {'background': tint(color, NODE_HIGHLIGHT_TINT), 'border': color}}


def style_network(net):
    '''
    Applies the look of the network: the label font, the node fills and the edge
    colours.
    '''
    net.options.nodes = {
        'font': {'size': LABEL_FONT_SIZE, 'color': LABEL_FONT_COLOR,
                 'face': LABEL_FONT_FACE,
                 # the halo that lifts the label off the edges
                 'strokeWidth': 4, 'strokeColor': NETWORK_BACKGROUND},
        'borderWidth': 2, 'borderWidthSelected': 3}
    net.options.edges = {
        # vis.js spreads the scores over 1-15 px by default, which turns the best
        # interactions into bars
        'scaling': {'min': 1, 'max': 6},
        # pyvis asks for dynamic curves, which hang a support node off every edge for the
        # physics
        'smooth': {'enabled': True, 'type': 'continuous', 'roundness': 0.15},
        # the edges are thin, so hovering or picking one thickens it
        'hoverWidth': 2, 'selectionWidth': 3}
    net.options.interaction = {'hover': True, 'selectConnectedEdges': False,
                               'tooltipDelay': 120}
    for node in net.nodes:
        node['color'] = node_color(node.get('color', '#8899aa'))
    for edge in net.edges:
        # as an object: vis.js reads a plain string as the colour whatever happens to the
        # edge
        edge['color'] = {'color': EDGE_COLOR,
                         'highlight': EDGE_ACCENT_COLOR,
                         'hover': EDGE_ACCENT_COLOR}


def generate_node_labels(df, annotations):
    '''
    Draws the descriptive protein name on the node instead of the short name STRING
    prefers.
    '''
    labels = {}
    for prefix in ['source', 'target']:
        for protein, name in df[[prefix, f'{prefix}_name']].drop_duplicates(subset=prefix).values:
            description = annotations.get(protein, '')
            if not description or description.lower().startswith('uncharacterized'):
                labels[protein] = str(name)
            else:
                labels[protein] = '\n'.join(textwrap.wrap(description, width=LABEL_WRAP_WIDTH))

    return labels


@st.cache_data(show_spinner=False)
def get_surface_calls(data_dir, host_taxids):
    '''
    What DeepLoc 2 called each protein and how sure it was of that call, keyed by STRING
    id.
    '''
    localisations = web_utils.load_deeploc_localisations(data_dir)
    if localisations.empty:
        return pd.DataFrame(columns=['surface'] + list(web_utils.DEEPLOC_SCORES.values()))

    hosts = localisations['protein'].str.split('.').str[0].isin(host_taxids)
    surface = pd.Series(index=localisations.index, dtype=object)
    surface[hosts] = web_utils.classify_localisation(
        localisations[hosts], web_utils.HOST_CLASSES, web_utils.SEVERAL)
    surface[~hosts] = web_utils.classify_localisation(
        localisations[~hosts], web_utils.SURFACE_CLASSES, web_utils.BOTH_SURFACE)
    scores = {column: localisations[column].values
              for column in web_utils.DEEPLOC_SCORES.values() if column in localisations}

    return pd.DataFrame({'surface': surface.values, **scores},
                        index=localisations['protein'].values)


def generate_node_titles(df, annotations, surface_calls=None):
    '''
    Builds the hover text of each node: the short name, the descriptive protein name,
    where DeepLoc puts the protein and the identifiers needed to look it up elsewhere.
    '''
    calls = {} if surface_calls is None or surface_calls.empty else surface_calls.to_dict('index')
    titles = {}
    # each side is read on the classes its own filter was made of
    for prefix, taxid_col, classes in [('source', 'taxid1_label', web_utils.SURFACE_CLASSES),
                                       ('target', 'taxid2_label', web_utils.HOST_CLASSES)]:
        cols = [prefix, f'{prefix}_name', f'{prefix}_uniprot', taxid_col]
        for protein, name, uniprot, species in df[cols].drop_duplicates(subset=prefix).values:
            lines = [str(name)]
            description = annotations.get(protein)
            if description and description != name:
                lines.append(description)
            lines.append(str(species))
            call = calls.get(protein)
            if call:
                # every class the protein was called for, each at its own probability
                called = [(c, call[web_utils.DEEPLOC_SCORES[c]]) for c in classes
                          if web_utils.DEEPLOC_SCORES[c] in call
                          and call[web_utils.DEEPLOC_SCORES[c]] > web_utils.DEEPLOC_CUTOFFS[c]]
                if called:
                    lines.append('DeepLoc: '
                                 + ', '.join(f'{c} (p={p:.2f})' for c, p in called))
                else:
                    lines.append(f"DeepLoc: {call['surface']}")
            lines.append(f'STRING: {protein}')
            if pd.notna(uniprot):
                lines.append(f'UniProt: {uniprot}')
            titles[protein] = '\n'.join(lines)

    return titles


def generate_interactions_table(df, score, annotations):
    '''Builds the table of predicted interactions above the chosen score.'''
    table = df[df['weight'] >= score]
    if table.empty:
        return table[TABLE_COLUMNS[:1]].rename(columns=TABLE_COLUMN_NAMES)

    tissues = table.groupby(['source', 'target'])['Tissue'].apply(
        lambda t: ', '.join(sorted(t.dropna().unique()))).rename('Tissues').reset_index()
    table = pd.merge(table.drop_duplicates(subset=['source', 'target']), tissues,
                     on=['source', 'target'])
    # the descriptive protein name, of which the predictions only keep the short version
    table = table.assign(
        source_full_name=table['source'].map(annotations).fillna(''),
        target_full_name=table['target'].map(annotations).fillna(''))

    columns = list(TABLE_COLUMNS)
    if table['taxid2_label'].nunique() > 1:
        columns.insert(columns.index('target_name'), 'taxid2_label')
    # where DeepLoc puts each side; only when the data directory has the localisations. The
    # parasite call is not a criterion
    if 'source_surface' in table.columns:
        columns.insert(columns.index('target_name'), 'source_surface')
    if 'target_surface' in table.columns:
        columns.insert(columns.index('Tissues'), 'target_surface')

    return (table[columns].sort_values(by='weight', ascending=False)
                          .rename(columns=TABLE_COLUMN_NAMES))


def search_table(table, search):
    '''
    Keeps the rows of the interactions table holding the searched text in any of their
    columns, ignoring case.
    '''
    search = search.strip()
    if not search or table.empty:
        return table

    matches = table.astype(str).apply(
        lambda column: column.str.contains(search, case=False, regex=False))

    return table[matches.any(axis=1)]


def name_the_selection(table, parasite, host):
    '''The parasite and the host as two columns at the front of the table.'''
    if table.empty:
        return table

    # the Rodent view already carries each row's host species
    host_column = 'Host group' if 'Host species' in table.columns else 'Host species'
    named = table.assign(**{'Parasite species': parasite, host_column: host})
    front = ['Parasite species', host_column]

    return named[front + [c for c in named.columns if c not in front]]


def active_filters(score, tissues, cell_types, surface):
    '''The filters narrowing what the page shows, as (kind, name) pairs.'''
    active = [('score', f'confidence \u2265 {score:g}')]
    for kind, selected in (('tissue', tissues), ('cell type', cell_types),
                           ('localisation', surface)):
        active.extend((kind, str(name)) for name in selected)

    return active


def show_active_filters(filters):
    '''Writes the filters above a section, so it says what it is showing.'''
    badges = ' '.join(f':{FILTER_BADGES.get(kind, "gray")}-badge[{name}]'
                      for kind, name in filters)
    st.markdown(f'Showing {badges}')


def download_name(parasite, filters, suffix):
    '''
    Names a download after the filters it was taken under, rather than after the
    parasite alone: an exported table of one tissue is otherwise indistinguishable from
    the whole network of that parasite once it is on disk.
    '''
    parts = [parasite]
    for kind in ('score', 'tissue', 'cell type', 'localisation'):
        names = [name for active_kind, name in filters if active_kind == kind]
        if not names:
            continue
        if kind == 'score':
            parts.append(names[0].split()[-1])
        elif len(names) <= FILENAME_MAX_NAMES:
            parts.extend(names)
        else:
            parts.append(f'{len(names)} {kind}s')

    slug = '_'.join(re.sub(r'[^A-Za-z0-9.]+', '-', part).strip('-') for part in parts)

    return f'{slug}_{suffix}'


def generate_tissue_filters(df):
    options = df['Tissue'].unique().tolist()

    return options

def generate_cell_type_filters(df, score):
    '''
    The cell types on offer, the one holding the most host proteins first, and how many
    each holds.
    '''
    annotated = df[(df['weight'] >= score) & df['Cell type'].notna()]
    if annotated.empty:
        return pd.Series(dtype=int)

    return (web_utils.keep_expressed_cell_types(annotated)
                     .groupby('Cell type')['target_name'].nunique()
                     .sort_values(ascending=False, kind='stable'))

def generate_surface_filters(df, surface_calls, niche):
    '''
    The DeepLoc classes offered as tickboxes: the surface of the host cell and the space
    around it, and for a parasite with an intracellular stage the cytosol and the
    nucleus as well, those being the classes its niche let the filter keep a host
    protein for.
    '''
    if surface_calls is None or surface_calls.empty or 'target' not in df.columns:
        return []
    called = surface_calls.reindex(df['target'].dropna().unique())

    return [c for c in web_utils.niche_classes(niche)
            if web_utils.DEEPLOC_SCORES[c] in called
            and (called[web_utils.DEEPLOC_SCORES[c]] > web_utils.DEEPLOC_CUTOFFS[c]).any()]


def cell_type_marks(df, score):
    '''
    The host proteins of the network against the cell types where they exceed 1 nTPM
    (web_utils.keep_expressed_cell_types), which is what both figures of the cell type
    section are drawn from, and the tissues those cell types are grouped into.
    '''
    annotated = df[(df['weight'] >= score) & df['Cell type'].notna()]
    if annotated.empty:
        return None, None

    marks = web_utils.keep_expressed_cell_types(annotated).copy()
    # the row a protein reaches its maximum in is always kept, so the darkest of a row is
    # 100%
    marks['share'] = marks['nTPM'] / marks.groupby(['target', 'Tissue'])['nTPM'].transform('max')
    marks['column'] = marks['Tissue'] + MATRIX_SEPARATOR + marks['Cell type']
    if (marks['target_name'].nunique() < MATRIX_MIN_PROTEINS
            or marks['column'].nunique() < 2):
        return None, None

    held = marks.groupby(['Tissue', 'column'])['target_name'].nunique()
    blocks = [(tissue, list(held[tissue].sort_values(ascending=False, kind='stable').index))
              for tissue in (marks.groupby('Tissue')['target_name'].nunique()
                                  .sort_values(ascending=False, kind='stable').index)]

    return marks, blocks


def label_tissue_blocks(figure, blocks):
    '''
    Names each tissue above the columns of its cell types and parts one block from the
    next with a rule, rather than repeating the tissue under every column of it.
    '''
    columns = [column for _, block in blocks for column in block]
    start = 0
    for tissue, block in blocks:
        if start:
            figure.add_vline(x=start - 0.5, line_width=1, line_color=GO_AXIS_COLOR)
        figure.add_annotation(x=start + (len(block) - 1) / 2, y=1.0, yref='paper',
                              yanchor='bottom', text=tissue, showarrow=False,
                              font=dict(color=LABEL_FONT_COLOR, size=13))
        start += len(block)

    figure.update_xaxes(tickmode='array', tickvals=columns, tickangle=-60, title=None,
                        showgrid=False,
                        ticktext=[column.split(MATRIX_SEPARATOR)[1] for column in columns])

    return figure


def style_cell_type_figure(figure, blocks):
    '''
    The two cell type figures share their axes, their blocks and the room above the plot
    the names of those blocks need, which style_go_axes leaves only the margin of a
    plain figure for.
    '''
    figure = style_go_axes(label_tissue_blocks(figure, blocks))
    figure.update_layout(margin=dict(l=10, r=10, t=34, b=40))

    return figure


def generate_cell_type_bars(marks, blocks):
    '''
    The predicted interactions of the network counted per cell type, in the blocks of
    the tissue the cell types belong to: where the parasite is predicted to meet the
    host most often, read in one glance.
    '''
    columns = [column for _, block in blocks for column in block]
    counted = (marks.drop_duplicates(['column', 'source', 'target_name'])
                    .groupby(['column', 'Tissue', 'Cell type'], observed=True).size()
                    .rename('interactions').reset_index())

    figure = px.bar(counted, x='column', y='interactions',
                    category_orders={'column': columns},
                    custom_data=['Cell type', 'Tissue'])
    figure.update_traces(marker_color=EDGE_ACCENT_COLOR,
                         hovertemplate='<b>%{customdata[0]}</b> (%{customdata[1]})<br>'
                                       'Predicted interactions: %{y}<extra></extra>')
    figure.update_yaxes(title='predicted interactions')
    web_utils.count_ticks(figure, counted['interactions'].max(), axis='y')
    figure.update_layout(height=BARS_HEIGHT)

    return style_cell_type_figure(figure, blocks)


def generate_cell_type_matrix(marks, blocks):
    '''
    Where inside the tissue the host proteins of the network sit: a mark wherever a
    protein exceeds 1 nTPM in a cell type, the cell types along the bottom in blocks of
    the tissue they belong to and the host proteins up the side.
    '''
    columns = [column for _, block in blocks for column in block]
    drawn = marks.drop_duplicates(['target_name', 'column'])
    proteins = list(drawn.groupby('target_name')['column'].nunique()
                         .sort_values(ascending=False, kind='stable').index)

    figure = px.scatter(drawn, x='column', y='target_name', color='share',
                        color_continuous_scale=GO_SEQUENTIAL,
                        range_color=(0, 1),
                        # plotly express flips category_orders on a y axis
                        category_orders={'column': columns, 'target_name': proteins},
                        custom_data=['Tissue', 'Cell type', 'nTPM', 'share'])
    figure.update_traces(marker=dict(size=MATRIX_MARK_SIZE, symbol='square',
                                     line=dict(color=NETWORK_BACKGROUND, width=1)),
                         hovertemplate='<b>%{y}</b><br>%{customdata[1]} (%{customdata[0]})'
                                       '<br>nTPM: %{customdata[2]:.1f}<br>'
                                       'Share of its peak in the tissue: '
                                       '%{customdata[3]:.0%}<extra></extra>')
    figure.update_yaxes(title=None, ticksuffix='  ')
    figure.update_coloraxes(colorbar=dict(title='Share of<br>tissue peak', tickformat='.0%',
                                          thickness=12, outlinewidth=0, len=0.6))
    # one row per host protein
    figure.update_layout(height=max(320, 22 * len(proteins) + 200))

    return style_cell_type_figure(figure, blocks)


@st.cache_data(max_entries=3, ttl=1800)
def get_enrichment(pred_df, data_dir, side, config_file):
    '''
    The processes over-represented among one side of the network, against the proteins
    of that side's species the network could have been drawn from.
    '''
    column, taxid_column = (('target', 'taxid2') if side == HOST else ('source', 'taxid1'))
    species = [int(s) for s in pred_df[taxid_column].unique()]
    # fastparquet prunes row groups only, so the exact selection is still applied afterwards
    go_df = utils.read_parquet_file(input_file=f'{data_dir}/gos.parquet', filters=[('taxid', 'in', species)])
    go_df = go_df[go_df['taxid'].isin(species)]
    # the background pool is read from the species of the view and the niche of its
    # parasite; several parasites at once are left on the union
    niche = None
    if side == HOST:
        niches = {web_utils.parasite_niche(utils.read_config(config_file), taxid)
                  for taxid in pred_df['taxid1'].unique()}
        niche = niches.pop() if len(niches) == 1 else None
    pool = web_utils.filtered_pool(data_dir, tuple(str(s) for s in species), niche=niche)
    # a directory carrying neither table is left on the proteome rather than emptied
    if pool:
        go_df = go_df[go_df['#string_protein_id'].isin(pool)]
    enrichment = utils.calculate_enrichment(set(pred_df[column]), go_df)
    # A is the number of proteins of the side annotated to the term
    enrichment = enrichment.rename(columns={'A': 'n_proteins'})

    return enrichment


def prepare_enrichment_view(enrichment_df):
    '''
    The columns the enrichment figures are drawn from: a plottable odds ratio (Fisher's
    exact test returns an infinite ratio for a term whose proteins are all in the
    network) and significance as -log10(FDR), which is what makes the significant terms
    spread out instead of piling up against zero.
    '''
    view = enrichment_df.copy()
    odds = pd.to_numeric(view['odds_ratio'], errors='coerce')
    finite = odds[np.isfinite(odds)]
    # an infinite ratio is drawn at the largest finite one and said so in the hover
    cap = finite.max() if not finite.empty else 1.0
    view['capped'] = ~np.isfinite(odds)
    view['odds_ratio_plot'] = odds.where(np.isfinite(odds), cap).clip(lower=np.nextafter(0, 1))
    fdr = pd.to_numeric(view['fdr_bh'], errors='coerce')
    smallest = fdr[fdr > 0].min()
    floor = smallest if pd.notna(smallest) else GO_MIN_FDR
    view['significance'] = -np.log10(fdr.clip(lower=max(floor, GO_MIN_FDR)))

    return view


def wrap_term(term):
    '''A GO term name broken over as many lines as it takes to fit beside an axis.'''
    return '<br>'.join(textwrap.wrap(str(term), width=GO_LABEL_WRAP_WIDTH)) or str(term)


def style_go_axes(fig):
    '''
    The grid and axes of the enrichment figures are a background to read the marks
    against, not lines to be read themselves.
    '''
    fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)',
                      font_color=LABEL_FONT_COLOR,
                      margin=dict(l=10, r=10, t=10, b=40), hovermode='closest')
    fig.update_xaxes(gridcolor=GO_GRID_COLOR, zeroline=False, linecolor=GO_AXIS_COLOR,
                     ticks='outside', tickcolor=GO_AXIS_COLOR, automargin=True)
    fig.update_yaxes(gridcolor=GO_GRID_COLOR, zeroline=False, linecolor=GO_AXIS_COLOR,
                     automargin=True)

    return fig


def get_enrichment_dotplot(enrichment_df, top_n=GO_TOP_N):
    '''
    The enriched processes as a ranked dot plot: the most significant terms, named on
    the axis, placed by odds ratio, sized by how many proteins of the network carry them
    and shaded by significance.
    '''
    view = prepare_enrichment_view(enrichment_df)
    view = view.nsmallest(top_n, 'fdr_bh').sort_values('odds_ratio_plot')
    view['label'] = view['go_term'].map(wrap_term)
    view['odds_ratio_text'] = np.where(view['capped'], 'infinite (all proteins in the network)',
                                       view['odds_ratio_plot'].map('{:.1f}'.format))

    fig = px.scatter(view, x='odds_ratio_plot', y='label', size='n_proteins',
                     color='significance', color_continuous_scale=GO_SEQUENTIAL,
                     size_max=26,
                     custom_data=['go_term', 'odds_ratio_text', 'n_proteins', 'fdr_bh'])
    fig.update_traces(marker=dict(line=dict(color=NETWORK_BACKGROUND, width=1.5)),
                      hovertemplate='<b>%{customdata[0]}</b><br>'
                                    'Odds ratio: %{customdata[1]}<br>'
                                    'Proteins in the network: %{customdata[2]}<br>'
                                    'FDR: %{customdata[3]:.2e}<extra></extra>')
    fig.update_xaxes(type='log', title='Odds ratio')
    fig.update_yaxes(title=None, showgrid=True, ticksuffix='  ')
    fig.update_coloraxes(colorbar=dict(title='-log<sub>10</sub> FDR', thickness=12,
                                       outlinewidth=0, len=0.6))
    # one row per process
    fig.update_layout(height=max(260, 90 + 34 * len(view)))

    return style_go_axes(fig)


def get_enrichment_volcano(enrichment_df, fdr, selected_terms=None, label_n=5):
    '''
    Every tested process, significant or not: effect size against significance, with the
    terms that pass the chosen FDR picked out of the ones that do not.
    '''
    selected_terms = set(selected_terms or [])
    view = prepare_enrichment_view(enrichment_df)
    view['log2_odds'] = np.log2(view['odds_ratio_plot'])
    view['significant'] = view['fdr_bh'] <= fdr
    view['odds_ratio_text'] = np.where(view['capped'], 'infinite (all proteins in the network)',
                                       view['odds_ratio_plot'].map('{:.1f}'.format))

    fig = pgo.Figure()
    # scaled by area from a reference shared by the three groups
    sizeref = 2.0 * view['n_proteins'].max() / (17.0 ** 2)
    groups = [('Not significant', view[~view['significant']], MUTED_COLOR),
              ('Significant', view[view['significant'] & ~view['go_term'].isin(selected_terms)],
               EDGE_ACCENT_COLOR),
              ('Selected', view[view['go_term'].isin(selected_terms)], HIGHLIGHT_COLOR)]
    for name, group, color in groups:
        if group.empty:
            continue
        fig.add_trace(pgo.Scatter(
            x=group['log2_odds'], y=group['significance'], mode='markers', name=name,
            marker=dict(color=color, size=group['n_proteins'], sizemode='area',
                        sizeref=sizeref, sizemin=4,
                        line=dict(color=NETWORK_BACKGROUND, width=1.2)),
            customdata=np.stack([group['go_term'], group['odds_ratio_text'],
                                 group['n_proteins'], group['fdr_bh']], axis=-1),
            hovertemplate='<b>%{customdata[0]}</b><br>'
                          'Odds ratio: %{customdata[1]}<br>'
                          'Proteins in the network: %{customdata[2]}<br>'
                          'FDR: %{customdata[3]:.2e}<extra></extra>'))

    # the most significant terms are named on the plot
    labelled = view[view['go_term'].isin(selected_terms)] if selected_terms else \
        view[view['significant']].nsmallest(label_n, 'fdr_bh')
    for i, (_, row) in enumerate(labelled.head(label_n).iterrows()):
        # staggered left and right so neighbouring labels do not overlap
        side = 1 if i % 2 == 0 else -1
        fig.add_annotation(x=row['log2_odds'], y=row['significance'],
                           text=textwrap.shorten(row['go_term'], width=34, placeholder='…'),
                           showarrow=True, arrowhead=0, arrowwidth=1,
                           arrowcolor=GO_AXIS_COLOR, ax=24 * side, ay=-16 - 12 * (i % 3),
                           xanchor='left' if side > 0 else 'right',
                           font=dict(size=11, color=LABEL_FONT_COLOR), align='left')

    fig.add_hline(y=-np.log10(fdr), line_dash='dot', line_color=GO_AXIS_COLOR,
                  annotation_text=f'FDR {fdr}', annotation_position='top left',
                  annotation_font=dict(size=11, color=GO_AXIS_COLOR))
    fig.update_xaxes(title='log<sub>2</sub> odds ratio')
    fig.update_yaxes(title='-log<sub>10</sub> FDR', rangemode='nonnegative')
    fig.update_layout(height=450, legend=dict(orientation='h', yanchor='bottom', y=1.02,
                                              x=0, title=None))

    return style_go_axes(fig)


@st.cache_data(show_spinner=False)
def load_ontology_parents(data_dir):
    '''child GO id -> its parent GO ids, the shape the ontology is walked upward in.'''
    ontology = load_ontology(data_dir)

    return ontology.groupby('child')['parent'].apply(list).to_dict()


def nearest_enriched_ancestors(terms, parents_of):
    '''
    For each enriched term, the closest term above it in the ontology that is also
    enriched, so the processes can be nested in each other.
    '''
    enriched = set(terms)
    ancestors = {}
    for term in terms:
        seen = {term}
        frontier = [p for p in parents_of.get(term, []) if p != term]
        found = ''
        while frontier:
            hits = sorted({p for p in frontier if p in enriched})
            if hits:
                found = hits[0]
                break
            next_frontier = []
            for node in frontier:
                if node in seen:
                    continue
                seen.add(node)
                next_frontier.extend(parents_of.get(node, []))
            frontier = [n for n in next_frontier if n not in seen]
        ancestors[term] = found

    # a cycle in the ontology would make the treemap undrawable, so the loop is cut at the
    # top
    for term in terms:
        walked = {term}
        node = ancestors[term]
        while node:
            if node in walked:
                ancestors[term] = ''
                break
            walked.add(node)
            node = ancestors.get(node, '')

    return ancestors


def picked(enrichment_df, selected_rows):
    '''
    Which rows of the enrichment the rows picked out of the grid are: by process, and by
    side when the grid shows both, since one process can be enriched on either.
    '''
    if SIDE_COLUMN in selected_rows:
        picks = set(zip(selected_rows[SIDE_COLUMN], selected_rows[GO_TERM_COLUMN]))
        return pd.Series(list(zip(enrichment_df['side'], enrichment_df['go_term'])),
                         index=enrichment_df.index).isin(picks)

    return enrichment_df['go_term'].isin(set(selected_rows[GO_TERM_COLUMN]))


def get_enrichment_summary(enrichment_df, parents_of):
    '''
    The enriched processes nested in the ontology: each term sits inside the closest
    enriched term above it, and the area of a block is the number of proteins of the
    network annotated to it.
    '''
    view = prepare_enrichment_view(enrichment_df).drop_duplicates(subset='go_id')
    terms = view['go_id'].tolist()
    labels = view['go_term'].tolist()
    ancestors = nearest_enriched_ancestors(terms, parents_of)
    odds_text = np.where(view['capped'], 'infinite',
                         view['odds_ratio_plot'].map('{:.1f}'.format))
    fdr_text = view['fdr_bh'].map('{:.2e}'.format)
    significance = view['significance']
    # the darker half of the ramp is too dark to write on in ink
    span = significance.max() - significance.min()
    midpoint = significance.min() + span / 2 if span > 0 else np.inf
    text_colors = np.where(significance > midpoint, '#ffffff', LABEL_FONT_COLOR)

    # the root gets the palest step of the ramp, and the scale is pinned to the terms
    fig = pgo.Figure(pgo.Treemap(
        # held by GO id and named by term, so two processes cannot collide
        ids=[GO_TREEMAP_ROOT_ID] + terms,
        # unwrapped: a block drops its name when it does not fit
        labels=[GO_TREEMAP_ROOT_LABEL] + labels,
        parents=[''] + [ancestors[t] or GO_TREEMAP_ROOT_ID for t in terms],
        values=[0] + view['n_proteins'].tolist(),
        # a parent term keeps its own proteins on top of those nested in it
        branchvalues='remainder',
        marker=dict(colors=[significance.min()] + significance.tolist(),
                    colorscale=GO_SEQUENTIAL,
                    cmin=significance.min(), cmax=significance.max(),
                    line=dict(color=NETWORK_BACKGROUND, width=2),
                    colorbar=dict(title='-log<sub>10</sub> FDR', thickness=12,
                                  outlinewidth=0, len=0.6)),
        customdata=[[GO_TREEMAP_ROOT_LABEL, '--', '--', '--']] +
                   np.stack([view['go_term'], odds_text,
                             view['n_proteins'].astype(str), fdr_text], axis=-1).tolist(),
        hovertemplate='<b>%{customdata[0]}</b><br>'
                      'Odds ratio: %{customdata[1]}<br>'
                      'Proteins in the network: %{customdata[2]}<br>'
                      'FDR: %{customdata[3]}<extra></extra>',
        insidetextfont=dict(color=[LABEL_FONT_COLOR] + text_colors.tolist()),
        # three levels of processes at a time (the root is the fourth); the rest by clicking
        # in
        maxdepth=4,
        pathbar=dict(visible=True, side='top', thickness=22),
        tiling=dict(pad=2)))
    fig.update_layout(height=700, margin=dict(l=0, r=0, t=10, b=10),
                      font_color=LABEL_FONT_COLOR, paper_bgcolor='rgba(0,0,0,0)',
                      # a name that does not fit is left out rather than shrunk
                      uniformtext=dict(minsize=9, mode='hide'))

    return fig

def generate_graph(df, score, annotations=None, surface_calls=None):
    if df.empty:
        return nx.Graph()

    G = nx.from_pandas_edgelist(df, 'source', 'target', 'weight')
    colors = dict(df[['source', 'source_color']].drop_duplicates().values)
    colors.update(dict(df[['target', 'target_color']].drop_duplicates().values))
    nx.set_node_attributes(G, colors, 'color')
    labels = dict(df[['source', 'source_name']].drop_duplicates().values)
    labels.update(dict(df[['target', 'target_name']].drop_duplicates().values))
    nx.set_node_attributes(G, labels, 'label')
    shapes = dict(df[['source', 'source_shape']].drop_duplicates().values)
    shapes.update(dict(df[['target', 'target_shape']].drop_duplicates().values))
    nx.set_node_attributes(G, shapes, 'shape')
    if annotations is not None:
        nx.set_node_attributes(G, generate_node_labels(df, annotations), 'label')
        nx.set_node_attributes(G, generate_node_titles(df, annotations, surface_calls), 'title')
    centrality = nx.betweenness_centrality(G, weight='weight')
    max_centrality = max(centrality.values(), default=0)
    sizes = {}
    for k,v in centrality.items():
        value = v*60/max_centrality if max_centrality > 0 else 20
        if value < 20:
            value = 20
        sizes[k] =  value
    nx.set_node_attributes(G, sizes, 'size')
    

    rm_edges = [(n1, n2) for n1,n2,w in G.edges.data('weight') if w < score]
    G.remove_edges_from(rm_edges)
    G.remove_nodes_from(list(nx.isolates(G)))

    widths = {}
    for n1,n2,w in df[['source', 'target', 'weight']].values:
        value = w*0.5/0.9
        if value < 0.05:
            value = 0.05
        widths[(n1, n2)] = value
    nx.set_edge_attributes(G, widths, 'value')
    nx.set_edge_attributes(G, EDGE_COLOR, 'color')


    return G


def annotate_edges(edges, df, annotations):
    '''
    Writes the two proteins of an interaction onto the edge that draws it, so that a
    click on the edge carries everything needed to show their structures.
    '''
    # the species travels with the edge, so the dialog can name it beside each protein
    cols = ['source', 'source_name', 'source_uniprot', 'taxid1_label',
            'target', 'target_name', 'target_uniprot', 'taxid2_label', 'weight',
            'group1', 'group2', 'experimental_evidence_score', 'databases_evidence_score']
    pairs = {}
    for row in df[cols].drop_duplicates(subset=['source', 'target']).itertuples(index=False):
        pairs[frozenset((row.source, row.target))] = {
            'parasite': str(row.source_name),
            # NaN does not survive the trip to the browser
            'parasite_uniprot': None if pd.isna(row.source_uniprot) else str(row.source_uniprot),
            'parasite_full': annotations.get(row.source, ''),
            'parasite_species': str(row.taxid1_label),
            'host': str(row.target_name),
            'host_uniprot': None if pd.isna(row.target_uniprot) else str(row.target_uniprot),
            'host_full': annotations.get(row.target, ''),
            'host_species': str(row.taxid2_label),
            'weight': float(row.weight),
            'parasite_group': str(row.group1),
            'host_group': str(row.group2),
            # the two scores the weight is the mean of, kept apart
            'experimental': float(row.experimental_evidence_score),
            'databases': float(row.databases_evidence_score)}

    annotated = []
    for edge in edges:
        edge = dict(edge)
        interaction = pairs.get(frozenset((edge['from'], edge['to'])))
        if interaction is not None:
            edge.update(interaction)
            # nothing else says the edges can be clicked
            edge['title'] = (f"{interaction['parasite']} -- {interaction['host']}\n"
                             f"confidence {interaction['weight']:.3f}\n"
                             'click to see the AlphaFold models')
        annotated.append(edge)

    return annotated


def network_options(net):
    '''The vis.js options pyvis built, which style_network has already filled in.'''
    return json.loads(net.get_network_data()[5])


@st.cache_data(show_spinner=False)
def get_structures(query_proteins):
    '''The AlphaFold model of each of the two proteins.'''
    return strv.get_alphafold_structure(query_proteins=query_proteins)


def show_structure(pdb_file):
    xyzview = strv.generate_mol_structure(pdb_file=pdb_file, height=VIEWER_HEIGHT)
    # same as stmol.showmol, which embeds through the deprecated st.components.v1.html
    st.iframe(xyzview._make_html(), height=VIEWER_HEIGHT + 20)


@st.dialog('AlphaFold models of the interacting proteins', width='large')
def show_structures_dialog(edge):
    '''
    Shows the AlphaFold model of each of the two proteins of the interaction that was
    clicked in the network, over the orthology groups the interaction was transferred
    from.
    '''
    st.markdown(f"**{edge['parasite']}** ({edge['parasite_species']}) &ndash; "
                f"**{edge['host']}** ({edge['host_species']}) &nbsp;·&nbsp; "
                f"interaction confidence score {edge['weight']:.2f}", unsafe_allow_html=True)
    st.caption(f"Transferred from the STRING link between orthology groups "
               f"**{edge['parasite_group']}** (parasite) and **{edge['host_group']}** "
               f"(host), scored on experimental evidence {edge['experimental']:.2f} and "
               f"database evidence {edge['databases']:.2f}.")
    st.markdown(strv.plddt_legend(), unsafe_allow_html=True)

    # one entry per protein, so the species stays attached even if the two proteins share a
    # name
    proteins = [(edge['parasite'], edge['parasite_uniprot'],
                 edge['parasite_full'], edge['parasite_species']),
                (edge['host'], edge['host_uniprot'],
                 edge['host_full'], edge['host_species'])]
    query_proteins = {edge['parasite']: edge['parasite_uniprot'],
                      edge['host']: edge['host_uniprot']}
    with st.spinner('Fetching the AlphaFold models...'):
        structures = get_structures(query_proteins)

    cols = st.columns(2)
    for i, (protein, uniprot, full_name, species) in enumerate(proteins):
        pdb_file, url, website, reason = structures[protein]
        with cols[i % len(cols)]:
            st.markdown(f'''<h4>{protein} ({species})</h4>''', unsafe_allow_html=True)
            subtitle = f'{full_name} · {uniprot}' if full_name else str(uniprot)
            st.caption(subtitle)
            if pdb_file is not None:
                show_structure(pdb_file=pdb_file)
                bcol1, bcol2 = st.columns(2)
                with bcol1:
                    st.link_button('PDB file', url)
                with bcol2:
                    st.link_button('AlphaFold EBI', website)
            else:
                st.markdown('''<h5>No AlphaFold model</h5>''', unsafe_allow_html=True)
                st.caption(reason)


st.caption('The predicted interactions of one host and one parasite as a network, with '
           'the AlphaFold model of the two proteins of an interaction and the biological '
           'processes over-represented among the host proteins of the network.')


col1, col2 = st.columns([1, 1], gap='large')

with col2:

    # the host carries over from whichever page it was last chosen on
    selected_host, selected_taxids = web_utils.host_selector(
        config, load_predictions(data_dir), 'Select a host', include_rodents=True)

    if selected_host == web_utils.NO_HOST:
        st.text('Choose 1 host to explore the predicted host-parasite interactions')
        selected_parasite = "<select>"
    else:
        # only parasites that infect the selected host
        parasite_list = ['<select>'] + get_parasite_list(data_dir, selected_taxids)
        # switching host can leave a parasite selected that the new host does not have
        if st.session_state.get('net_par') not in parasite_list:
            st.session_state.pop('net_par', None)
        selected_parasite = st.selectbox('Select a parasite to visualize the predicted PPI', parasite_list, key="net_par")

    if selected_parasite != "<select>":
        df_select = get_parasite_tissues(data_dir, selected_parasite, selected_taxids)
        df_select = web_utils.filter_tissues(config, df_select)
        # carried on the predictions so the filter, the table and the network read the same
        # call
        surface_calls = get_surface_calls(data_dir, tuple(str(t) for t in selected_taxids))
        if not surface_calls.empty:
            df_select = df_select.assign(
                source_surface=df_select['source'].map(surface_calls['surface']),
                target_surface=df_select['target'].map(surface_calls['surface']))
        score = st.slider('Confidence score', 0.35, 0.9, 0.35)

        tissues_options = generate_tissue_filters(df_select)
        if len(tissues_options) > 0:
            # a click on the body figure drawn below is read here, on the run after it,
            # before the filter is created
            body_figure.apply_organ_click(config, selected_taxids, tissues_options,
                                          TISSUE_FILTER_KEY)
            selected_tissues = st.multiselect('Select tissues to filter the predicted PPI',
                                              tissues_options, key=TISSUE_FILTER_KEY)
            if len(selected_tissues) > 0:
                df_select = df_select[df_select['Tissue'].isin(selected_tissues)]

        # cell types are offered on their own; picking tissues first narrows what is offered
        cell_type_counts = generate_cell_type_filters(df_select, score)
        if len(cell_type_counts) > 0:
            def cell_type_label(cell_type):
                proteins = cell_type_counts[cell_type]

                return f'{cell_type} ({proteins} protein{"" if proteins == 1 else "s"})'

            selected_cell_types = st.multiselect(
                'Select cell types to filter the predicted PPI', list(cell_type_counts.index),
                format_func=cell_type_label,
                help='A cell type is offered with the number of host proteins expressed '
                     'above 1 nTPM in it. Cell-type annotation is available for human (HPA) '
                     'and pig (Pig Cell Atlas); a host protein with no cell type is left out '
                     'once a cell type is chosen.')
            if len(selected_cell_types) > 0:
                expressed = web_utils.keep_expressed_cell_types(
                    df_select[df_select['Cell type'].notna()])
                df_select = expressed[expressed['Cell type'].isin(selected_cell_types)]

        # the tickboxes are drawn above the network; Streamlit hands their state over before
        # they are drawn again
        surface_options = generate_surface_filters(
            df_select, surface_calls,
            web_utils.get_niches(config).get(selected_parasite, web_utils.UNKNOWN_NICHE))
        ticked = [c for c in surface_options
                  if st.session_state.get(SURFACE_FILTER_KEYS[c])]
        # ticking every class leaves every host protein in, as ticking none does
        selected_surface = ticked if 0 < len(ticked) < len(surface_options) else []
        if selected_surface:
            # read from the probabilities so a protein called for several answers to each
            over = pd.concat([surface_calls[web_utils.DEEPLOC_SCORES[c]]
                              > web_utils.DEEPLOC_CUTOFFS[c]
                              for c in selected_surface], axis=1).any(axis=1)
            df_select = df_select[df_select['target'].map(over).fillna(False)]

        annotations = web_utils.load_protein_annotations(data_dir)
        for host_taxid in selected_taxids:
            host_df = df_select[df_select['taxid2'].astype(str) == host_taxid]
            host_df = host_df.assign(
                target_color=config['hosts'][int(host_taxid)]['color'])
            G = generate_graph(host_df, score, annotations, surface_calls)
            net = Network(height=f'{NETWORK_HEIGHT}px', width="100%",
                          bgcolor=NETWORK_BACKGROUND, font_color=LABEL_FONT_COLOR)
            net.from_nx(G)
            net.repulsion(node_distance=420, central_gravity=0.33,
                          spring_length=110, spring_strength=0.10,
                          damping=0.95)
            style_network(net)
            networks.append((host_taxid, config['hosts'][int(host_taxid)]['label'],
                             host_df, G, net))
        
        
        
# resolved after the column that draws the filters, read by every section below
page_filters = active_filters(score, selected_tissues, selected_cell_types,
                              selected_surface) if df_select is not None else []

# drawn after the column that holds the selectors, where the predictions are filtered
with col1:
    if df_select is not None:
        body_figure.show_body_figure(config, data_dir, df_select[df_select['weight'] >= score],
                                     selected_taxids, selected_tissues, clickable=True)


def network_legend(parasite_label, parasite_color, host_label, host_color):
    '''
    The key to the network: a node is read by its shape, which says whether the protein
    is the parasite's or the host's, and by its colour, which says which species it
    belongs to.
    '''
    def entry(shape, color, label):
        fill = tint(color, NODE_FILL_TINT)
        mark = (f'<polygon points="11,2 20,11 11,20 2,11"'
                if shape == 'diamond' else f'<circle cx="11" cy="11" r="8.5"')

        return (f'<span style="display:inline-flex;align-items:center;gap:0.45em;">'
                f'<svg width="22" height="22" viewBox="0 0 22 22">'
                f'{mark} fill="{fill}" stroke="{color}" stroke-width="2"/></svg>'
                f'<span>{label}</span></span>')

    return (f'<div style="display:flex;flex-wrap:wrap;align-items:center;gap:1.6em;'
            f'margin:0.2rem 0 0.6rem;font-size:0.8rem;color:{LABEL_FONT_COLOR};">'
            f'{entry("diamond", parasite_color, parasite_label)}'
            f'{entry("dot", host_color, host_label)}</div>')


def render_network_panel(host_taxid, host_label, host_df, G, net):
    if net is not None:
        # the species name only distinguishes the panels on the rodent page
        if len(networks) > 1:
            st.subheader(host_label)
        st.text(f"Nodes: {len(G.nodes())}  Edges: {len(G.edges())}")
        if host_df.empty:
            st.info(f'No predicted interactions for this parasite in {host_label}.')
            return
        st.markdown(network_legend(host_df['taxid1_label'].iloc[0],
                                   host_df['source_color'].iloc[0],
                                   host_label, host_df['target_color'].iloc[0]),
                    unsafe_allow_html=True)
        filename = f'{selected_parasite}_{host_taxid}_network'
        html_data = ""
        # saved for the download button; the page draws the network through the component so
        # clicks reach Python
        net.save_graph(f'{path}/{filename}.html')
        utils.export_graph(G, filename=f'{filename}.graphml',
                           format='graphml', output_dir=path)
        utils.export_graph(G, filename=f'{filename}.json',
                           format='cytoscape', output_dir=path)
        with open(f'{path}/{filename}.html','r',encoding='utf-8') as HtmlFile:
            html_data = HtmlFile.read()
        nodes, edges = net.get_network_data()[:2]
        selected_edge = ppi_network(
            nodes=nodes,
            edges=annotate_edges(edges, host_df, web_utils.load_protein_annotations(data_dir)),
            options=network_options(net),
            height=NETWORK_HEIGHT,
            key=f'network_{selected_parasite}_{host_taxid}')
        net = None

        # closing the dialog reruns the page with the component still holding the edge, so
        # the click is remembered; the nonce lets the same edge open again
        shown_edge_key = f'shown_edge_{host_taxid}'
        if selected_edge is not None and selected_edge['nonce'] != st.session_state.get(shown_edge_key):
            st.session_state[shown_edge_key] = selected_edge['nonce']
            show_structures_dialog(selected_edge['edge'])
        with st.container():
            c1, c2, c3 = st.columns(3)

            with c1:
                st.download_button(
                    label="Download Network as Html",
                    data=html_data,
                    file_name=download_name(selected_parasite, page_filters,
                                            f'{host_taxid}_network.html'),
                    mime='text/html',
                )
            with c2:
                st.download_button(
                    label="Download Network as GraphML",
                    data=open(f'{path}/{filename}.graphml','r',encoding='utf-8'),
                    file_name=download_name(selected_parasite, page_filters,
                                            f'{host_taxid}_network.graphml'),
                    mime='text/plain',
                )
            with c3:
                st.download_button(
                    label="Download Network as Cytoscape",
                    data=open(f'{path}/{filename}.json','r',encoding='utf-8'),
                    file_name=download_name(selected_parasite, page_filters,
                                            f'{host_taxid}_network.json'),
                    mime='text/plain',
                )


if networks:
    st.header('Network of host-parasite PPIs')
    show_active_filters(page_filters)
    st.caption('Predicted interactions between parasite and host proteins above the '
               'selected confidence score. Nodes are proteins, diamonds parasite and '
               'circles host, coloured by organism and sized by centrality in the '
               'network. Edge width is the confidence score. Hover a node for its full '
               'name and identifiers; click an edge for the AlphaFold model of both '
               'proteins.')
    if selected_host == 'Rodent':
        st.caption('Rat and Mouse networks are shown separately; their host nodes use '
                   'their species-specific colors.')
    if len(surface_options) > 0:
        st.caption('The host proteins are in the predictions because DeepLoc called them '
                   'where DeepLoc puts them: on the membrane of the host cell, '
                   'extracellular -- in the matrix and the fluid around it -- or, for a '
                   'parasite with an intracellular stage, in the cytosol or the nucleus. '
                   'Ticking a class leaves the interactions that can take place there, and '
                   'drops any parasite protein left with nothing to bind; a protein DeepLoc '
                   'places in several classes stays whichever of them is ticked.')
        # one column more than there are boxes so they are not spread across the page
        boxes = st.columns(len(surface_options) + 2)
        for box, surface in zip(boxes, surface_options):
            with box:
                st.checkbox(surface, key=SURFACE_FILTER_KEYS[surface])
    columns = st.columns(len(networks))
    for column, network in zip(columns, networks):
        with column:
            render_network_panel(*network)


with st.container():
    if df_select is not None:
        marks, blocks = cell_type_marks(df_select, score)
        if marks is not None:
            st.header('Cell types expressing the host proteins')
            show_active_filters(page_filters)
            st.caption('A host protein counts towards a cell type when its expression is '
                       'above 1 nTPM. The columns of both tabs are those cell types, grouped '
                       'into the tissues the parasite infects, and a cell type is written '
                       'under its block alone, since the same kind of cell is annotated '
                       'separately in each tissue. Cell-type annotation is available for '
                       'human (HPA) and pig (Pig Cell Atlas).')
            # the same columns counted, then opened up per protein
            per_cell_type_tab, per_protein_tab = st.tabs(['Per cell type', 'Per protein'])
            with per_cell_type_tab:
                st.caption('Predicted interactions per cell type. An interaction is counted '
                           'in every cell type where its host protein exceeds 1 nTPM, so the '
                           'bars overlap and are not a partition of the network.')
                st.plotly_chart(generate_cell_type_bars(marks, blocks), width='stretch')
            with per_protein_tab:
                st.caption('A mark wherever a host protein exceeds 1 nTPM in a cell type, '
                           'shaded by the share of its expression in that tissue the cell '
                           'type carries. A row of one mark is a protein the parasite meets '
                           'in a single kind of cell; a full row one it meets throughout the '
                           'tissue.')
                st.plotly_chart(generate_cell_type_matrix(marks, blocks), width='stretch')

with st.container():
    if df_select is not None:
        st.header("Table of host-parasite PPIs")
        show_active_filters(page_filters)
        table = generate_interactions_table(df_select, score,
                                            web_utils.load_protein_annotations(data_dir))
        st.caption('One row per predicted interaction, with the tissues the host protein is '
                   'expressed in and its DeepLoc class.')
        search = st.text_input(
            "Search the table",
            key='table_search',
            placeholder="Protein name, identifier, tissue ...",
            help="Keeps the rows holding the text typed, looked for in any column.")
        table = search_table(table, search)
        if search.strip() and table.empty:
            st.info(f"No interaction holds '{search}'.")
        gb = GridOptionsBuilder.from_dataframe(table)
        # a page of a known size with the grid grown to fit it, rather than blank rows under
        # a short page
        gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=TABLE_PAGE_SIZE)
        gb.configure_grid_options(domLayout='autoHeight')
        gb.configure_side_bar()
        gridOptions = gb.build()
        grid_response = AgGrid(
                            table,
                            gridOptions=gridOptions,
                            data_return_mode='AS_INPUT',
                            fit_columns_on_grid_load=False,
                            enable_enterprise_modules=True
                        )
        st.download_button(
            label="Download Network Table",
            data=utils.convert_df(name_the_selection(table, selected_parasite,
                                                     selected_host)),
            file_name=download_name(selected_parasite, page_filters, 'network_table.tsv'),
            mime='text/csv',
        )

with st.container():
    if df_select is not None:
        st.header("Functional enrichment of the network (GO biological processes)")
        show_active_filters(page_filters)
        st.caption('Biological processes over-represented among the host proteins, the '
                   'parasite proteins, or both, each side tested on its own. '
                   'A side is tested against the proteins of its own species the '
                   'pipeline had to work with -- the ones its filters passed, the host '
                   'proteins on expression and localisation and the parasite proteins on '
                   'being secreted -- and not against the whole proteome, which would '
                   'return those filters themselves as a result. The two sides are tested '
                   'apart: they are annotated to a different depth and were selected on '
                   "different grounds. One-sided Fisher's exact test, corrected across "
                   'terms with Benjamini-Hochberg; a process is tested when the background '
                   'gives it at least 11 proteins, no more than a quarter of them, and at '
                   'least two are in the network.')
        side = st.radio('Proteins to test', (HOST, PARASITE, BOTH), horizontal=True,
                        help='The host proteins the parasite is predicted to reach, the '
                             'parasite proteins reaching them, or the two tested apart and '
                             'shown together.')
        sides = [HOST, PARASITE] if side == BOTH else [side]
        # the sides are tested one at a time and carried together with the side on each row
        enrichment = pd.concat([get_enrichment(df_select[df_select['weight'] >= score],
                                               data_dir, s,
                                               web_utils.get_config_file()).assign(side=s)
                                for s in sides], ignore_index=True)
        if not enrichment.empty:
            fdr = st.radio('False discovery rate', (0.01, 0.05, 0.1), index=1,
                           horizontal=True,
                           help='The Benjamini-Hochberg corrected significance a process '
                                'has to reach to be counted as enriched.')
            # only the grid and the file it hands out are renamed for reading
            enrichment_view = enrichment[enrichment['fdr_bh'] <= fdr]
            # the side is a column of the grid only when both are in it
            column_names = ({'side': SIDE_COLUMN} if side == BOTH else {}) | ENRICHMENT_COLUMN_NAMES
            enrichment_table = enrichment_view[list(column_names)].rename(columns=column_names)
            # counted per side when both are shown
            by_side = enrichment_view['side'].value_counts()
            per_side = (' (' + ', '.join(f'{by_side.get(s, 0)} among the {s.lower()}'
                                         for s in sides) + ')') if side == BOTH else ''
            st.caption(f'{len(enrichment_table)} processes pass an FDR of {fdr}{per_side}. '
                       'Select rows to highlight GO terms in the network.')
            gb = GridOptionsBuilder.from_dataframe(enrichment_table)
            gb.configure_pagination(paginationAutoPageSize=False,
                                    paginationPageSize=TABLE_PAGE_SIZE)
            gb.configure_grid_options(domLayout='autoHeight')
            gb.configure_side_bar()
            gb.configure_selection('multiple', use_checkbox=True, groupSelectsChildren="Group checkbox select children")
            gridOptions = gb.build()
            grid_response = AgGrid(
                                enrichment_table,
                                gridOptions=gridOptions,
                                data_return_mode='AS_INPUT',
                                fit_columns_on_grid_load=False,
                                enable_enterprise_modules=True
                            )
            selected_rows = grid_response['selected_rows']
            st.download_button(
                label="Download Enrichment Table",
                data=utils.convert_df(enrichment_table),
                file_name=f'{selected_parasite}_network_enrichment_table.tsv',
                mime='text/csv',
            )
        else:
            st.subheader("No GO terms were found enriched")

with st.container():
    if enrichment_view is not None and enrichment_view.empty:
        st.info(f"No biological process passes an FDR of {fdr}. Loosen the correction above "
                "to see the processes that are enriched less strongly.")
    elif enrichment_view is not None:
        enrichment_viz = enrichment_view
        if selected_rows is not None and len(selected_rows) > 0:
            selected_terms = selected_rows[GO_TERM_COLUMN].values.tolist()
            enrichment_viz = enrichment_viz[picked(enrichment_viz, selected_rows)]

        st.subheader("Enriched biological processes")
        ranked_tab, volcano_tab = st.tabs(["Ranked processes", "All tested processes"])
        # one figure per side shown, each named when there are two
        with ranked_tab:
            st.caption(f'The {GO_TOP_N} most significantly over-represented biological '
                       f'processes among the {"each side" if side == BOTH else side.lower()} '
                       'of the network, positioned by odds ratio, sized by the number of '
                       'them annotated to each process and shaded by significance. Select '
                       'rows in the table above to restrict the figure.')
            for s in sides:
                part = enrichment_viz[enrichment_viz['side'] == s]
                if side == BOTH:
                    st.markdown(f'**{s}**')
                if part.empty:
                    st.caption(f'No process passes an FDR of {fdr} among the {s.lower()}.')
                    continue
                st.plotly_chart(get_enrichment_dotplot(part), width='stretch',
                                key=f'dotplot_{s}')
        with volcano_tab:
            st.caption('All processes tested against the network: effect size on the x axis '
                       'and significance on the y axis, with the processes passing the '
                       'selected FDR highlighted.')
            for column, s in zip(st.columns(len(sides)), sides):
                part = enrichment[enrichment['side'] == s]
                with column:
                    if side == BOTH:
                        st.markdown(f'**{s}**')
                    if part.empty:
                        st.caption(f'No process could be tested among the {s.lower()}.')
                        continue
                    picked_terms = enrichment_viz[enrichment_viz['side'] == s]['go_term'] \
                        .tolist() if selected_terms else []
                    st.plotly_chart(get_enrichment_volcano(part, fdr, picked_terms),
                                    width='stretch', key=f'volcano_{s}')

        with st.container():
            if len(selected_terms) > 0:
                if enrichment is not None:
                    highlighted_nodes = enrichment[picked(enrichment, selected_rows)]['nodes'].values
                    highlighted_nodes = utils.merge_list_of_lists([i.split(',') for i in highlighted_nodes])
                    highlight_color = {i: HIGHLIGHT_COLOR for i in highlighted_nodes}
                    G = generate_graph(df_select, score,
                                       web_utils.load_protein_annotations(data_dir),
                                       surface_calls)
                    nx.set_node_attributes(G, MUTED_COLOR, 'color')
                    nx.set_node_attributes(G, highlight_color, 'color')
                    net = Network(height="450px", width="100%",
                                  bgcolor=NETWORK_BACKGROUND, font_color=LABEL_FONT_COLOR)
                    net.from_nx(G)
                    G = None
                    style_network(net)
                    net.save_graph(f'{path}/{selected_parasite}2.html')
                    net = None
                    st.subheader("Nodes annotated to the selected biological processes")
                    st.caption('The network with the proteins annotated to the biological '
                               'processes selected above in pink and the remainder in grey.')
                    with open(f'{path}/{selected_parasite}2.html','r',encoding='utf-8') as HtmlFile:
                        html_data = HtmlFile.read()
                    st.iframe(html_data, height=500)
                    st.download_button(
                        label="Download Network as Html",
                        data=html_data,
                        file_name=f'{selected_parasite}_enrichment_network.html',
                        mime='text/html',
                    )
        
        st.subheader("Hierarchy of enriched biological processes")
        st.caption('Enriched processes arranged by the Gene Ontology hierarchy: each process '
                   'is nested within the closest enriched process above it, its area is the '
                   'number of network proteins annotated to it and its shade is its '
                   'significance. Click a block to zoom in.')
        # one hierarchy per side, named when both are shown
        for s in sides:
            part = enrichment_view[enrichment_view['side'] == s]
            if side == BOTH:
                st.markdown(f'**{s}**')
            if part.empty:
                st.caption(f'No process passes an FDR of {fdr} among the {s.lower()}.')
                continue
            st.plotly_chart(get_enrichment_summary(part, load_ontology_parents(data_dir)),
                            width='stretch', key=f'treemap_{s}')

st.markdown("---")
st.markdown("---")

# Footer
with st.container():
    web_utils.footer()