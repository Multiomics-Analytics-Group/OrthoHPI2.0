import sys, os
import json
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

#Initialize variables
df_select = None
net = None
selected_rows = []
selected_terms = []
enrichment_table = None
enrichment_view = None
enrichment = None
path = 'data/tmp'
# where the network is written for the download buttons; not in the repository, so it
# has to be created on a fresh checkout
os.makedirs(path, exist_ok=True)
# characters per line of a node label, so a long protein name wraps instead of
# stretching its node across the network
LABEL_WRAP_WIDTH = 22
LABEL_FONT_SIZE = 22  # vis.js defaults to 14, too small to read the protein names
LABEL_FONT_COLOR = '#2f3a46'
LABEL_FONT_FACE = "'Source Sans Pro', -apple-system, 'Segoe UI', sans-serif"
# the labels are drawn over the edges, so they are given a halo of the network
# background to sit in rather than being read through the lines
NETWORK_BACKGROUND = '#fbfcfd'
# how far toward white a node's species colour is washed out to fill it. A saturated
# fill of the colours the parasites are identified by leaves the label on top of it
# unreadable; a wash of the same colour still says which organism the protein is from
NODE_FILL_TINT = 0.74
NODE_HIGHLIGHT_TINT = 0.5
# edges are the background of the picture until one is pointed at: neutral and light
# enough to read the nodes through them, and the accent of the page when hovered
EDGE_COLOR = '#c7cfd9'
EDGE_ACCENT_COLOR = '#2b8cbe'
# the enrichment view below the table, where the proteins of the selected processes are
# picked out of a network that is otherwise pushed into the background
HIGHLIGHT_COLOR = '#e7298a'
MUTED_COLOR = '#c8ced6'
# the enrichment figures. Significance is a magnitude, so it is drawn on a single hue
# running light to dark rather than on a set of unrelated colours
GO_SEQUENTIAL = ['#bcdcec', '#7fc0dd', '#3f9fca', '#2b8cbe', '#12587d', '#08324a']
GO_AXIS_COLOR = '#8d97a3'
# the matrix of host proteins against the cell types they are concentrated in. Fewer
# proteins than this and the rows say nothing a glance at the table would not
MATRIX_MIN_PROTEINS = 3
MATRIX_MARK_SIZE = 11
# the bars are counted per cell type and do not grow with the host proteins, so the tab
# holding them keeps one height whatever the parasite
BARS_HEIGHT = 420
# what a column joins the tissue and the cell type with, since a cell type name is only
# unique inside its tissue. Not a character a name of either carries
MATRIX_SEPARATOR = ' | '
GO_GRID_COLOR = '#e6eaef'
# how many processes the ranked plot shows. Beyond this the term names stop being
# readable, and the table above is the place to see the rest
GO_TOP_N = 20
# characters per line of a GO term name on the y axis of the ranked plot
GO_LABEL_WRAP_WIDTH = 42
# an FDR that underflows to 0 would be an infinite -log10. It is read as the smallest
# FDR the rest of the table holds instead, so one term cannot stretch the colour scale
GO_MIN_FDR = 1e-300
# the block all the others are nested in. It is in the data rather than left to plotly,
# which paints a sector it has no colour value for in a hard grey (see the treemap below)
GO_TREEMAP_ROOT_ID = '__all_enriched_processes__'
GO_TREEMAP_ROOT_LABEL = 'All enriched processes'
# height of the network. The layout spaces the proteins widely so that their names can be
# read, and the view is zoomed to fit whatever it is given: a shorter canvas does not
# remove the space between the proteins, it only shrinks the whole network, names and all.
# The empty space around a small network is taken out by the zoom (see index.html), not
# by making the canvas smaller.
NETWORK_HEIGHT = 1000
# height of each of the two structure viewers in the dialog an interaction opens
VIEWER_HEIGHT = 400

# Columns of the interactions table. The colors and shapes only exist to draw the network,
# the taxids repeat their labels, and the parasite and the edge type are the same on every
# row. The selected host is one species, so its proteins are unambiguous throughout the
# network, table, and tissue views.
#
# Both STRING ids are carried even though the parasite one is the taxid and the UniProt
# accession joined and could be left to the reader to assemble: the host one cannot be --
# 9606.ENSP00000312671 shares nothing with Q8N3D4 -- and a table whose two halves are
# identified differently is a table that has to be explained before it can be used.
TABLE_COLUMNS = ['source_name', 'source_full_name', 'source', 'source_uniprot',
                 'target_name', 'target_full_name', 'target', 'target_uniprot',
                 'Tissues', 'experimental_evidence_score',
                 'databases_evidence_score', 'weight', 'group1', 'group2']

# what those columns are called once they are read by someone rather than by the network:
# `source` is the parasite and `target` the host throughout the pipeline, which is not
# something the table should ask its reader to know. The names are the ones the rest of the
# app uses, so a column here and a caption elsewhere say the same word for the same thing.
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

# and the columns of the enrichment table. GO_TERM_COLUMN is read back out of the grid
# when rows are picked, so the figures below highlight what was selected.
ENRICHMENT_COLUMN_NAMES = {
    'go_term': 'Biological process',
    'n_proteins': 'Proteins of the network',
    'odds_ratio': 'Odds ratio',
    'p_value': 'P value',
    'fdr_bh': 'FDR (BH)'}
GO_TERM_COLUMN = ENRICHMENT_COLUMN_NAMES['go_term']

# rows to a page of either table. The grid is grown to the page rather than the page fitted
# to a fixed height, so the last page of a short table is the only one drawn short
TABLE_PAGE_SIZE = 10

# the two sides of the network, tested separately for enrichment
HOST, PARASITE = 'Host proteins', 'Parasite proteins'

# what the host proteins of the network are tested against. The pipeline's filters are
# the universe the network was drawn from, and the second option narrows that universe
# the same way the predictions themselves are narrowed. The two answer different
# questions -- whether the parasite's targets are special among everything it could have
# reached anywhere, or among what it could have reached where it lives -- and the wider
# one is offered first because it is the same background for every parasite, so two
# parasites can be read against each other
BACKGROUND_FILTERS = 'Proteins that passed the filters'
BACKGROUND_TISSUES = 'Those in the tissues this parasite infects'

# Read dataset
config = utils.read_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()


# the predictions themselves are loaded by web_utils, so the three pages share one
# cached copy of them however the app is navigated
load_predictions = web_utils.load_predictions


@st.cache_data(show_spinner=False, max_entries=5)
def get_parasite_tissues(data_dir, parasite, host_taxids):
    '''
    Annotates the predictions of one parasite against the selected host with the
    tissues and cell types its host targets are expressed in. Selecting the parasite
    before merging keeps the full predictions x tissues table (~630 MB) out of the app.
    '''
    predictions = load_predictions(data_dir)
    tissues = utils.read_parquet_file(input_file=f'{data_dir}/tissues_cell_types.parquet')
    predictions = predictions[(predictions['taxid1_label'] == parasite)
                              & (predictions['taxid2'].isin(host_taxids))]

    return pd.merge(predictions, tissues.rename({'Gene': 'target'}, axis=1), on='target', how='left')


@st.cache_data(show_spinner=False)
def get_parasite_list(data_dir, host_taxids):
    '''Parasites with predictions against the selected host, so the parasite
    selector never offers one that yields an empty network.'''
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

    :param str color: '#rrggbb'
    :param float amount: 0 leaves the colour alone, 1 turns it white
    :return: the mixed colour, or the colour unchanged if it is not a hex triplet
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
    a group at a glance and keeps the protein name on top of it legible. Hovering and
    selecting deepen the fill rather than change the colour, so a node stays
    recognisable as its species while it is pointed at.

    :param str color: the species colour
    :return: a vis.js node colour dictionary
    '''
    return {'background': tint(color, NODE_FILL_TINT),
            'border': color,
            'highlight': {'background': tint(color, NODE_HIGHLIGHT_TINT), 'border': color},
            'hover': {'background': tint(color, NODE_HIGHLIGHT_TINT), 'border': color}}


def style_network(net):
    '''
    Applies the look of the network: the label font, the node fills and the edge
    colours. The font has to be set on the network rather than on each node, as pyvis
    overwrites a node's font with the font_color the network was built with. The colours
    are set on the nodes and edges pyvis built rather than on the networkx graph they
    came from, because a colour dictionary is not something GraphML can hold and the
    same graph is exported for download.

    :param net: pyvis Network to style
    '''
    net.options.nodes = {
        'font': {'size': LABEL_FONT_SIZE, 'color': LABEL_FONT_COLOR,
                 'face': LABEL_FONT_FACE,
                 # the halo that lifts the label off the edges running under it
                 'strokeWidth': 4, 'strokeColor': NETWORK_BACKGROUND},
        'borderWidth': 2, 'borderWidthSelected': 3}
    net.options.edges = {
        # vis.js spreads the confidence scores over its default 1-15 px, which turns the
        # best-supported interactions into bars. A narrower range still ranks them by
        # confidence while leaving the network something to be read through
        'scaling': {'min': 1, 'max': 6},
        # pyvis asks for dynamic curves, which hang an invisible support node off every
        # edge for the physics to solve: the same curve without the cost, and shallow
        # enough that two proteins still read as joined by a line
        'smooth': {'enabled': True, 'type': 'continuous', 'roundness': 0.15},
        # the edges are thin, which is close to unclickable, so pointing at one or
        # picking it thickens it as well as colouring it
        'hoverWidth': 2, 'selectionWidth': 3}
    net.options.interaction = {'hover': True, 'selectConnectedEdges': False,
                               'tooltipDelay': 120}
    for node in net.nodes:
        node['color'] = node_color(node.get('color', '#8899aa'))
    for edge in net.edges:
        # as an object rather than pyvis' plain string: vis.js reads a string as "this
        # colour whatever happens to the edge", which loses the hover and the selection
        edge['color'] = {'color': EDGE_COLOR,
                         'highlight': EDGE_ACCENT_COLOR,
                         'hover': EDGE_ACCENT_COLOR}


def generate_node_labels(df, annotations):
    '''
    Draws the descriptive protein name on the node instead of the short name STRING
    prefers. Proteins STRING has no name for are left with their short name -- an
    "Uncharacterized protein" label identifies nothing, and most parasite proteins
    would end up sharing it. The name is wrapped so long ones do not stretch the node
    across the network.

    :param df: predictions dataframe of the selected parasite
    :param dict annotations: STRING id --> descriptive protein name
    :return: {STRING id: label}
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
def get_surface_calls(data_dir):
    '''
    What DeepLoc 2 called each protein and how sure it was of that call, keyed by STRING
    id. Both sides of the interactions are in it: every protein of the predictions went
    through a localisation filter to get here, the host proteins for being surface-exposed
    and the parasite proteins for being secreted, and this is what those filters read.

    The score kept beside the class is the probability of that class -- the one the call
    was made on -- and not both probabilities, since it is the one a tooltip has room for.

    :param str data_dir: directory holding deeploc_localisations.parquet
    :return: dataframe of surface and score, indexed by STRING id; empty without the file
    '''
    localisations = web_utils.load_deeploc_localisations(data_dir)
    if localisations.empty:
        return pd.DataFrame(columns=['surface', 'score'])

    surface = web_utils.classify_surface(localisations)

    return pd.DataFrame({
        'surface': surface.values,
        'score': np.where(surface == web_utils.EXTRACELLULAR, localisations['extracellular'],
                          localisations['cell_membrane'])},
        index=localisations['protein'].values)


def generate_node_titles(df, annotations, surface_calls=None):
    '''
    Builds the hover text of each node: the short name, the descriptive protein name,
    where DeepLoc puts the protein and the identifiers needed to look it up elsewhere.
    The label only shows one of the two names, so the tooltip keeps both.

    :param df: predictions dataframe of the selected parasite
    :param dict annotations: STRING id --> descriptive protein name
    :param surface_calls: DeepLoc calls, as get_surface_calls returns them, or None for a
                          data directory with no localisations to show
    :return: {STRING id: hover text}
    '''
    calls = {} if surface_calls is None or surface_calls.empty else surface_calls.to_dict('index')
    titles = {}
    for prefix, taxid_col in [('source', 'taxid1_label'), ('target', 'taxid2_label')]:
        cols = [prefix, f'{prefix}_name', f'{prefix}_uniprot', taxid_col]
        for protein, name, uniprot, species in df[cols].drop_duplicates(subset=prefix).values:
            lines = [str(name)]
            description = annotations.get(protein)
            if description and description != name:
                lines.append(description)
            lines.append(str(species))
            call = calls.get(protein)
            if call:
                lines.append(f"DeepLoc: {call['surface']} (p={call['score']:.2f})")
            lines.append(f'STRING: {protein}')
            if pd.notna(uniprot):
                lines.append(f'UniProt: {uniprot}')
            titles[protein] = '\n'.join(lines)

    return titles


def generate_interactions_table(df, score, annotations):
    '''
    Builds the table of predicted interactions above the chosen score. The predictions
    are merged with the tissue and cell type annotation to filter them, which repeats
    every interaction once per tissue and single-cell cluster the host protein was
    found in -- 39 interactions come out as 820 rows. The tissues are gathered back
    into one cell so the table holds one row per interaction.

    :param df: predictions dataframe of the selected parasite
    :param float score: minimum confidence score
    :param dict annotations: STRING id --> descriptive protein name
    :return: the table to show and download
    '''
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
    # where DeepLoc puts each side. The host call is why the host protein is in the network
    # at all; the parasite call is not a criterion -- the secretome filter is what selected
    # those -- and is here to be read rather than to explain the selection, which is why a
    # parasite protein can be called neither class and still be in the table. Only when the
    # data directory has the localisations, since the snapshot ones do not
    if 'source_surface' in table.columns:
        columns.insert(columns.index('target_name'), 'source_surface')
    if 'target_surface' in table.columns:
        columns.insert(columns.index('Tissues'), 'target_surface')

    return (table[columns].sort_values(by='weight', ascending=False)
                          .rename(columns=TABLE_COLUMN_NAMES))


def search_table(table, search):
    '''
    Keeps the rows of the interactions table holding the searched text in any of their
    columns, ignoring case. The rows shown and the rows handed out by the download
    button are the same ones this way.

    :param table: the interactions table
    :param str search: text typed in the search box, empty for no filtering
    :return: the rows that match
    '''
    search = search.strip()
    if not search or table.empty:
        return table

    matches = table.astype(str).apply(
        lambda column: column.str.contains(search, case=False, regex=False))

    return table[matches.any(axis=1)]


def name_the_selection(table, parasite, host):
    '''
    The parasite and the host as two columns at the front of the table.

    For the table on the page they would be the same value on every row -- one parasite and
    one host group are selected before it is built -- and the host is only named there when
    the group covers two species. The downloaded file is the other way round: it is read
    away from the selectors that made it, and outlives the parasite name in its filename.

    :param table: the interactions table, as generate_interactions_table returns it
    :param str parasite: the selected parasite
    :param str host: the selected host species
    :return: the table with the two columns in front, or unchanged if it is empty
    '''
    if table.empty:
        return table

    # The table may already contain a host-species column from an upstream export; keep it
    # and label this selected-species column as the host group in that case.
    host_column = 'Host group' if 'Host species' in table.columns else 'Host species'
    named = table.assign(**{'Parasite species': parasite, host_column: host})
    front = ['Parasite species', host_column]

    return named[front + [c for c in named.columns if c not in front]]


def generate_tissue_filters(df):
    options = df['Tissue'].unique().tolist()

    return options

def generate_cell_type_filters(df, score):
    '''
    The cell types on offer, the one holding the most host proteins first, and how many
    each holds. The count is what makes one cell type worth picking over another, and it
    is read off the interactions above the confidence slider so that an option never
    promises proteins the network is not drawing.

    A cell type holds the host proteins concentrated in it rather than those merely
    detected there (web_utils.keep_peak_cell_types, which says why), and the filter that
    reads these options keeps the same rows.

    :param dataframe df: the predictions of the parasite, annotated with tissues
    :param float score: confidence the network is drawn from
    :return: series of host proteins per cell type, its index the options themselves
    '''
    annotated = df[(df['weight'] >= score) & df['Cell type'].notna()]
    if annotated.empty:
        return pd.Series(dtype=int)

    return (web_utils.keep_peak_cell_types(annotated)
                     .groupby('Cell type')['target_name'].nunique()
                     .sort_values(ascending=False, kind='stable'))

def generate_surface_filters(df):
    '''
    The DeepLoc classes the host proteins of this parasite fall in, in the order they are
    offered: the surface of the host cell first and the space around it second, which is
    the order the home page splits its columns in. A class no host protein of this parasite
    is in is left out, as an empty option filters to an empty network.
    '''
    if 'target_surface' not in df.columns:
        return []
    present = set(df['target_surface'].dropna())

    return [c for c in (web_utils.CELL_MEMBRANE, web_utils.EXTRACELLULAR) if c in present]


def cell_type_marks(df, score):
    '''
    The host proteins of the network against the cell types they are concentrated in
    (web_utils.keep_peak_cell_types), which is what both figures of the cell type section
    are drawn from, and the tissues those cell types are grouped into.

    A cell type name repeats across tissues -- smooth muscle cells are in the lung and in
    the intestine -- so a column is a (tissue, cell type) pair and only the cell type is
    written under it. The tissues come in blocks, the one holding the most host proteins
    first and its cell types in the same order inside it, as the tissues of the other pages
    are ordered.

    :param dataframe df: the predictions of the parasite, annotated with tissues
    :param float score: confidence the network is drawn from
    :return: the rows and the blocks as (tissue, its columns), or (None, None) where there
             is too little annotation to draw anything
    '''
    annotated = df[(df['weight'] >= score) & df['Cell type'].notna()]
    if annotated.empty:
        return None, None

    marks = web_utils.keep_peak_cell_types(annotated).copy()
    # the row a protein reaches its maximum in is always kept, so the share is read off
    # the marks themselves and the darkest of a row is 100%
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
    next with a rule, rather than repeating the tissue under every column of it. The two
    figures carry the same columns in the same order, so the blocks land in the same place
    on both and switching between them does not move the ground under the reader.

    A line down every column would fence the marks in and leave these rules
    indistinguishable from the rest, so the vertical grid goes as they arrive.

    :param figure: a figure whose x axis holds the columns of the blocks, in their order
    :param list blocks: the blocks as (tissue, its columns)
    :return: the figure
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
    '''The two cell type figures share their axes, their blocks and the room above the
    plot the names of those blocks need, which style_go_axes leaves only the margin of a
    plain figure for.'''
    figure = style_go_axes(label_tissue_blocks(figure, blocks))
    figure.update_layout(margin=dict(l=10, r=10, t=34, b=40))

    return figure


def generate_cell_type_bars(marks, blocks):
    '''
    The predicted interactions of the network counted per cell type, in the blocks of the
    tissue the cell types belong to: where the parasite is predicted to meet the host most
    often, read in one glance.

    An interaction is counted in every cell type its host protein is concentrated in, so
    the bars overlap and do not partition the network. They are the columns of the matrix
    beside them added up, which is the trade the two tabs offer: how many against which.

    :param dataframe marks: rows from cell_type_marks
    :param list blocks: the blocks as (tissue, its columns)
    :return: the figure
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
    Where inside the tissue the host proteins of the network sit: a mark wherever a protein
    is concentrated in a cell type, the cell types along the bottom in blocks of the tissue
    they belong to and the host proteins up the side.

    The network says which host proteins a parasite is predicted to reach and the body
    figure says in which organs. Neither can say whether a protein is met in one cell type
    of an organ or in all of them, which is what a row of this reads as: a row of a single
    mark is a protein the parasite meets in one kind of cell, a full row one it meets
    wherever it goes.

    The mark is shaded by the share of the expression the protein reaches anywhere in that
    tissue the cell type carries, so the darkest mark of a row is where the protein is most
    abundant.

    :param dataframe marks: rows from cell_type_marks
    :param list blocks: the blocks as (tissue, its columns)
    :return: the figure
    '''
    columns = [column for _, block in blocks for column in block]
    drawn = marks.drop_duplicates(['target_name', 'column'])
    proteins = list(drawn.groupby('target_name')['column'].nunique()
                         .sort_values(ascending=False, kind='stable').index)

    figure = px.scatter(drawn, x='column', y='target_name', color='share',
                        color_continuous_scale=GO_SEQUENTIAL,
                        range_color=(web_utils.PEAK_CELL_TYPE_FRACTION, 1),
                        # plotly express flips category_orders on a y axis, so the protein
                        # in the most cell types first puts it in the top row
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
    # one row is one host protein, as one row of the enrichment dot plot is one process
    figure.update_layout(height=max(320, 22 * len(proteins) + 200))

    return style_cell_type_figure(figure, blocks)


@st.cache_data(max_entries=3, ttl=1800)
def get_enrichment(pred_df, data_dir, side, background, config_file):
    '''
    The processes over-represented among one side of the network, against the proteins of
    that side's species the network could have been drawn from.

    The two sides are tested apart rather than pooled. They are annotated to a different
    depth -- 17,265 human proteins carry a term against 11,025 of S. mansoni -- and they
    were selected on different grounds, the host proteins for being surface or
    extracellular and the parasite proteins for being secreted. Pooled, the test is run
    against a null that is half the other organism, and whichever side has more proteins in
    the network decides what the section says.

    Either side is tested against the proteins of its species that came through the
    pipeline's filters, not against the whole proteome. Those filters are why a network
    holds the proteins it holds, so a proteome background reads them back as a result: at
    the lowest confidence the page offers, the whole S. mansoni network returns 383
    processes against the proteome, led by cell-substrate adhesion -- which is the surface
    call that let those host proteins in -- and 104 against the pool. Its parasite side
    goes from 49 processes to 28, tested against the 826 proteins the secretome filter
    passed rather than the 11,025 annotated ones. The difference is smaller on the
    networks a high confidence leaves, which run to a few proteins and can as easily gain
    terms as lose them: what changes there is which terms, not how many.

    :param pred_df: predictions of the network, above the chosen score
    :param str data_dir: directory holding gos.parquet
    :param str side: HOST or PARASITE
    :param str background: BACKGROUND_FILTERS or BACKGROUND_TISSUES, for the host side
    :param str config_file: configuration naming the tissues the parasite infects, part of
                            the cache key so a snapshot entrypoint tests against its own
    :return: enrichment dataframe, with A renamed to n_proteins
    '''
    column, taxid_column = (('target', 'taxid2') if side == HOST else ('source', 'taxid1'))
    species = [int(s) for s in pred_df[taxid_column].unique()]
    # the filter is pushed down to the reader (fastparquet prunes row groups only,
    # so the exact selection is still applied afterwards)
    go_df = utils.read_parquet_file(input_file=f'{data_dir}/gos.parquet', filters=[('taxid', 'in', species)])
    go_df = go_df[go_df['taxid'].isin(species)]
    # The selected host is one species, so the background pool is read from its own taxid.
    pool = web_utils.filtered_pool(data_dir, tuple(str(s) for s in species))
    if side == HOST and background == BACKGROUND_TISSUES:
        parasite = pred_df['taxid1'].iloc[0]
        pool = pool & web_utils.infected_tissue_proteins(
            data_dir, utils.read_config(config_file), parasite)
    # a data directory carrying neither table gives no pool, and is left on the proteome
    # rather than emptied -- the same fallback the snapshot directories need elsewhere
    if pool:
        go_df = go_df[go_df['#string_protein_id'].isin(pool)]
    enrichment = utils.calculate_enrichment(set(pred_df[column]), go_df)
    # A is the number of proteins of the side annotated to the term, which is what every
    # figure below sizes its marks by
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
    # an infinite ratio is drawn at the largest one that can be drawn, and said so in
    # the hover, rather than dropped: those terms are the strongest results
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
    '''The grid and axes of the enrichment figures are a background to read the marks
    against, not lines to be read themselves.'''
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
    the axis, placed by odds ratio, sized by how many proteins of the network carry
    them and shaded by significance.
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
    # one row of the plot is one process, so the canvas grows with the number of them
    # instead of squeezing the names together
    fig.update_layout(height=max(260, 90 + 34 * len(view)))

    return style_go_axes(fig)


def get_enrichment_volcano(enrichment_df, fdr, selected_terms=None, label_n=5):
    '''
    Every tested process, significant or not: effect size against significance, with
    the terms that pass the chosen FDR picked out of the ones that do not.
    '''
    selected_terms = set(selected_terms or [])
    view = prepare_enrichment_view(enrichment_df)
    view['log2_odds'] = np.log2(view['odds_ratio_plot'])
    view['significant'] = view['fdr_bh'] <= fdr
    view['odds_ratio_text'] = np.where(view['capped'], 'infinite (all proteins in the network)',
                                       view['odds_ratio_plot'].map('{:.1f}'.format))

    fig = pgo.Figure()
    # the marks are scaled by area from a reference shared by the three groups, so a dot
    # means the same number of proteins wherever it sits
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

    # the few most significant terms are named on the plot, so the shape of the cloud
    # can be read without pointing at it
    labelled = view[view['go_term'].isin(selected_terms)] if selected_terms else \
        view[view['significant']].nsmallest(label_n, 'fdr_bh')
    for i, (_, row) in enumerate(labelled.head(label_n).iterrows()):
        # the labels are staggered left and right so that the names of two terms of
        # nearly the same significance do not land on top of each other
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
    enriched, so the processes can be nested in each other. Terms with no enriched
    ancestor sit at the top level. The ontology is a graph rather than a tree -- a term
    can have several parents -- so it is walked breadth-first and the first enriched
    level reached wins.

    Terms are GO ids here rather than names: the ontology is keyed by id, and a name
    that does not match a node in it would leave the term looking like a root.
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

    # two terms can end up as each other's ancestor if the ontology holds a cycle;
    # a treemap cannot be drawn from one, so the loop is cut at the top
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
    # the darker half of the colour ramp is too dark to write the term name on in ink
    span = significance.max() - significance.min()
    midpoint = significance.min() + span / 2 if span > 0 else np.inf
    text_colors = np.where(significance > midpoint, '#ffffff', LABEL_FONT_COLOR)

    # the root is drawn behind every block, and plotly paints a sector it has no colour
    # value for in a hard grey -- the frame that used to sit around the whole treemap.
    # It is given the palest step of the ramp, and the scale is pinned to the terms so
    # that the root cannot shift it
    fig = pgo.Figure(pgo.Treemap(
        # a block is held by its GO id and named by its term, so two processes cannot
        # collide and the nesting follows the ontology rather than the wording
        ids=[GO_TREEMAP_ROOT_ID] + terms,
        # unwrapped: a block draws its name on one line and drops it when it does not
        # fit, which keeps the header of a group readable instead of hiding it
        labels=[GO_TREEMAP_ROOT_LABEL] + labels,
        parents=[''] + [ancestors[t] or GO_TREEMAP_ROOT_ID for t in terms],
        values=[0] + view['n_proteins'].tolist(),
        # a parent term keeps its own proteins on top of those of the terms nested in it
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
        # a network can be enriched for several hundred terms, which at full depth are
        # slivers too small to read. Three levels of processes are drawn at a time (the
        # root is the fourth) and the rest is reached by clicking into a block
        maxdepth=4,
        pathbar=dict(visible=True, side='top', thickness=22),
        tiling=dict(pad=2)))
    fig.update_layout(height=700, margin=dict(l=0, r=0, t=10, b=10),
                      font_color=LABEL_FONT_COLOR, paper_bgcolor='rgba(0,0,0,0)',
                      # a name that does not fit its block is left out rather than
                      # shrunk to an unreadable size
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
    # remove filtered edges from graph G
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
    click on the edge carries everything needed to show their structures. The graph is
    undirected, which leaves the orientation of an edge up to networkx -- some come out
    host to parasite -- so the pair is looked up unordered and written back with the
    parasite protein always first.

    What the prediction was transferred from travels with the edge too: the two orthology
    groups and the two evidence scores STRING gives their link. The edge is drawn between
    two proteins, but nothing about those two proteins was measured -- the pair exists
    because their groups interact -- so the groups are what a reader has to see to know
    what the edge stands on.

    :param list edges: vis.js edge dictionaries, as pyvis built them
    :param df: predictions dataframe of the selected parasite
    :param dict annotations: STRING id --> descriptive protein name
    :return: the same edges, each with the proteins of its interaction added
    '''
    # the species travels with the edge as well, so the dialog can name it beside each
    # protein and say which of the two models is the parasite's and which the host's
    cols = ['source', 'source_name', 'source_uniprot', 'taxid1_label',
            'target', 'target_name', 'target_uniprot', 'taxid2_label', 'weight',
            'group1', 'group2', 'experimental_evidence_score', 'databases_evidence_score']
    pairs = {}
    for row in df[cols].drop_duplicates(subset=['source', 'target']).itertuples(index=False):
        pairs[frozenset((row.source, row.target))] = {
            'parasite': str(row.source_name),
            # a protein UniProt has no accession for has no AlphaFold model either,
            # and NaN does not survive the trip to the browser
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
            # the two scores the weight is the mean of, kept apart: a link carried by
            # experiments and one carried by a database curation average out the same and
            # are not the same evidence
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
    '''
    The vis.js options pyvis built, which style_network has already filled in. The
    downloaded HTML is drawn from the same options, so the network that is saved looks
    like the one on the page.

    :param net: pyvis Network the options come from
    :return: options dictionary to hand to the network component
    '''
    return json.loads(net.get_network_data()[5])


@st.cache_data(show_spinner=False)
def get_structures(query_proteins):
    '''
    The AlphaFold model of each of the two proteins. Cached so that reopening a pair is
    immediate and a rerun of the page never downloads anything again.

    :param dict query_proteins: protein name --> UniProt accession
    :return: dict protein name --> (pdb file, pdb url, AlphaFold entry url)
    '''
    return strv.get_alphafold_structure(query_proteins=query_proteins)


def show_structure(pdb_file):
    xyzview = strv.generate_mol_structure(pdb_file=pdb_file, height=VIEWER_HEIGHT)
    # same as stmol.showmol, which still embeds through the deprecated st.components.v1.html.
    # the width is left to stretch so the viewer follows the column it sits in
    st.iframe(xyzview._make_html(), height=VIEWER_HEIGHT + 20)


@st.dialog('AlphaFold models of the interacting proteins', width='large')
def show_structures_dialog(edge):
    '''
    Shows the AlphaFold model of each of the two proteins of the interaction that was
    clicked in the network, over the orthology groups the interaction was transferred from.
    Opened over the network rather than placed under it so that the table and the
    enrichment plots below do not move on every click.

    The provenance is named as a pair of groups and not as a pair of proteins, which is
    what it is: the link comes from STRING's COG links, one row per pair of orthology
    groups, and no interaction between these two proteins was ever measured. Writing it
    the other way round -- an interaction between two named proteins somewhere in STRING
    -- would be inventing a record the transfer never had.

    :param dict edge: the clicked edge, as annotated by annotate_edges
    '''
    st.markdown(f"**{edge['parasite']}** ({edge['parasite_species']}) &ndash; "
                f"**{edge['host']}** ({edge['host_species']}) &nbsp;·&nbsp; "
                f"interaction confidence score {edge['weight']:.2f}", unsafe_allow_html=True)
    st.caption(f"Transferred from the STRING link between orthology groups "
               f"**{edge['parasite_group']}** (parasite) and **{edge['host_group']}** "
               f"(host), scored on experimental evidence {edge['experimental']:.2f} and "
               f"database evidence {edge['databases']:.2f}.")
    st.markdown(strv.plddt_legend(), unsafe_allow_html=True)

    # one entry per protein rather than a dict keyed by protein name, so the species stays
    # attached to the right panel even if the two proteins happen to share a name
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


# the body figure beside the selectors rather than a third of the row left empty
col1, col2 = st.columns([1, 1], gap='large')

with col2:

    # the host carries over from whichever page it was last chosen on
    selected_host, selected_taxids = web_utils.host_selector(
        config, load_predictions(data_dir), 'Select a host')

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

    # Set info message on initial site load
    if selected_parasite == "<select>":
        st.text('Choose 1 parasite to visualize the predicted PPI network')
    else:
        df_select = get_parasite_tissues(data_dir, selected_parasite, selected_taxids)
        df_select = web_utils.filter_tissues(config, df_select)
        # where DeepLoc puts each host protein, carried on the predictions so the filter
        # below, the table and the network all read the same call
        surface_calls = get_surface_calls(data_dir)
        if not surface_calls.empty:
            df_select = df_select.assign(
                source_surface=df_select['source'].map(surface_calls['surface']),
                target_surface=df_select['target'].map(surface_calls['surface']))
        score = st.slider('Confidence score', 0.4, 0.9, 0.4)

        selected_tissues = []
        tissues_options = generate_tissue_filters(df_select)
        if len(tissues_options) > 0:
            selected_tissues = st.multiselect('Select tissues to filter the predicted PPI', tissues_options)
            if len(selected_tissues) > 0:
                df_select = df_select[df_select['Tissue'].isin(selected_tissues)]

        # the cell types are offered on their own rather than only after a tissue has been
        # picked: a cell type belongs to one tissue anyway, so choosing one is choosing its
        # tissue as well, and someone after a cell type had to find its tissue first to be
        # allowed to ask for it. Picking tissues first narrows what is offered here, since
        # the options are read off whatever is left of the predictions
        cell_type_counts = generate_cell_type_filters(df_select, score)
        if len(cell_type_counts) > 0:
            def cell_type_label(cell_type):
                proteins = cell_type_counts[cell_type]

                return f'{cell_type} ({proteins} protein{"" if proteins == 1 else "s"})'

            selected_cell_types = st.multiselect(
                'Select cell types to filter the predicted PPI', list(cell_type_counts.index),
                format_func=cell_type_label,
                help='A cell type is offered with the number of host proteins concentrated '
                     'in it -- expressed there at half at least of what they reach anywhere '
                     'in the tissue, since nearly every protein of a tissue is detected in '
                     'nearly every one of its cell types. The single cell annotation is '
                     'human, so this is offered for human host proteins alone, and a host '
                     'protein with no cell type is left out once a cell type is chosen.')
            if len(selected_cell_types) > 0:
                peak = web_utils.keep_peak_cell_types(df_select[df_select['Cell type'].notna()])
                df_select = peak[peak['Cell type'].isin(selected_cell_types)]

        # localisation is independent of the tissue, so it filters beside the tissues
        # rather than inside them: a host protein is on the cell surface or in the space
        # around it wherever it is expressed
        surface_options = generate_surface_filters(df_select)
        if len(surface_options) > 0:
            selected_surface = st.multiselect(
                'Select where DeepLoc places the host proteins', surface_options,
                help='The host proteins are in the predictions because DeepLoc called them '
                     'surface-exposed: on the membrane of the host cell, or extracellular -- '
                     'in the matrix and the fluid around it. Filtering to one of the two '
                     'leaves the interactions that can take place there, and drops any '
                     'parasite protein left with nothing to bind.')
            if len(selected_surface) > 0:
                df_select = df_select[df_select['target_surface'].isin(selected_surface)]

        # Create networkx graph object from pandas dataframe
        G = generate_graph(df_select, score, web_utils.load_protein_annotations(data_dir),
                           surface_calls)
            
        st.text(f"Nodes: {len(G.nodes())}  Edges: {len(G.edges())}")

        # Initiate PyVis network object
        net = Network(height=f'{NETWORK_HEIGHT}px', width="100%",
                      bgcolor=NETWORK_BACKGROUND, font_color=LABEL_FONT_COLOR)
        # Take Networkx graph and translate it to a PyVis graph format
        net.from_nx(G)
        # Save other formats
        utils.export_graph(G, filename=f'{selected_parasite}.graphml',
                        format='graphml', output_dir=f'{path}')
        utils.export_graph(G, filename=f'{selected_parasite}.json',
                        format='cytoscape', output_dir=f'{path}')
        G = None

        # Generate network with specific layout settings. The repulsion solver rather
        # than force atlas: it spaces the proteins far enough apart that the names on
        # them can be read, which is worth more than the tighter, rounder shape the
        # force atlas layout draws.
        net.repulsion(node_distance=420, central_gravity=0.33,
                        spring_length=110, spring_strength=0.10,
                        damping=0.95)
        style_network(net)
        
        #net.show_buttons(filter_=['nodes'])
        
        
# drawn after the column that holds the selectors, which is where the predictions the
# figure counts are read and filtered
with col1:
    if df_select is not None:
        body_figure.show_body_figure(config, data_dir, df_select[df_select['weight'] >= score],
                                     selected_taxids, selected_tissues)


with st.container():
    if net is not None:
        st.header('Network of host-parasite PPIs')
        st.caption('Predicted interactions between parasite and host proteins above the '
                   'selected confidence score. Nodes are proteins, diamonds parasite and '
                   'circles host, coloured by organism and sized by centrality in the '
                   'network. Edge width is the confidence score. Hover a node for its full '
                   'name and identifiers; click an edge for the AlphaFold model of both '
                   'proteins.')
        html_data = ""
        # Save and read graph as HTML file, which is what the download button hands out.
        # The network on the page is drawn by the component instead: an embedded HTML
        # file has no way of telling Python which interaction was clicked.
        net.save_graph(f'{path}/{selected_parasite}.html')
        with open(f'{path}/{selected_parasite}.html','r',encoding='utf-8') as HtmlFile:
            html_data = HtmlFile.read()
        nodes, edges = net.get_network_data()[:2]
        # no widget key on purpose: Streamlit then identifies the component by its
        # arguments, so choosing another parasite or moving the score slider makes it a
        # different widget and the interaction selected in the previous network is
        # dropped rather than carried over to one that no longer contains it
        selected_edge = ppi_network(
            nodes=nodes,
            edges=annotate_edges(edges, df_select, web_utils.load_protein_annotations(data_dir)),
            options=network_options(net),
            height=NETWORK_HEIGHT)
        net = None

        # Closing the dialog reruns the page with the component still holding the edge
        # that opened it, so the click is remembered to keep it from opening again. The
        # nonce changes on every click, which is what makes clicking the same edge twice
        # open the dialog again.
        if selected_edge is not None and selected_edge['nonce'] != st.session_state.get('shown_edge'):
            st.session_state['shown_edge'] = selected_edge['nonce']
            show_structures_dialog(selected_edge['edge'])
        with st.container():
            c1, c2, c3 = st.columns(3)

            with c1:
                st.download_button(
                    label="Download Network as Html",
                    data=html_data,
                    file_name=f'{selected_parasite}_network.html',
                    mime='text/html',
                )
            with c2:
                st.download_button(
                    label="Download Network as GraphML",
                    data=open(f'{path}/{selected_parasite}.graphml','r',encoding='utf-8'),
                    file_name=f'{selected_parasite}_network.graphml',
                    mime='text/plain',
                )
            with c3:
                st.download_button(
                    label="Download Network as Cytoscape",
                    data=open(f'{path}/{selected_parasite}.json','r',encoding='utf-8'),
                    file_name=f'{selected_parasite}_network.json',
                    mime='text/plain',
                )


with st.container():
    if df_select is not None:
        marks, blocks = cell_type_marks(df_select, score)
        if marks is not None:
            st.header('Cell types the host proteins are concentrated in')
            st.caption('A host protein counts towards a cell type where it is concentrated: '
                       'expressed there at half at least of what it reaches anywhere in that '
                       'tissue. The columns of both tabs are those cell types, grouped into '
                       'the tissues the parasite infects, and a cell type is written under '
                       'its block alone, since the same kind of cell is annotated separately '
                       'in each tissue. The single cell annotation is human, so this is drawn '
                       'for human host proteins alone.')
            # the same columns twice, counted and then opened up: how many interactions a
            # cell type holds, and which host proteins they are
            per_cell_type_tab, per_protein_tab = st.tabs(['Per cell type', 'Per protein'])
            with per_cell_type_tab:
                st.caption('Predicted interactions per cell type. An interaction is counted '
                           'in every cell type its host protein is concentrated in, so the '
                           'bars overlap and are not a partition of the network.')
                st.plotly_chart(generate_cell_type_bars(marks, blocks), width='stretch')
            with per_protein_tab:
                st.caption('A mark wherever a host protein is concentrated in a cell type, '
                           'shaded by the share of its expression in that tissue the cell '
                           'type carries. A row of one mark is a protein the parasite meets '
                           'in a single kind of cell; a full row one it meets throughout the '
                           'tissue.')
                st.plotly_chart(generate_cell_type_matrix(marks, blocks), width='stretch')

with st.container():
    if df_select is not None:
        st.header("Table of host-parasite PPIs")
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
        # a page of a known size with the grid grown to fit it. Sizing the page to a fixed
        # height instead fits as many whole rows as it can and draws the remainder of the
        # height as blank rows under the last one
        gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=TABLE_PAGE_SIZE)
        gb.configure_grid_options(domLayout='autoHeight')
        gb.configure_side_bar() #Add a sidebar
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
            file_name=f'{selected_parasite}_network_table.tsv',
            mime='text/csv',
        )

with st.container():
    if df_select is not None:
        st.header("Functional enrichment of the network (GO biological processes)")
        st.caption('Biological processes over-represented among one side of the network. '
                   'Each side is tested against the proteins of its own species the '
                   'pipeline had to work with -- the ones its filters passed, the host '
                   'proteins on expression and localisation and the parasite proteins on '
                   'being secreted -- and not against the whole proteome, which would '
                   'return those filters themselves as a result. The two sides are tested '
                   'apart: they are annotated to a different depth and were selected on '
                   "different grounds. One-sided Fisher's exact test, corrected across "
                   'terms with Benjamini-Hochberg; a process is tested when the background '
                   'gives it at least 11 proteins, no more than a quarter of them, and at '
                   'least two are in the network.')
        side = st.radio('Proteins to test', (HOST, PARASITE), horizontal=True,
                        help='The host proteins the parasite is predicted to reach, or the '
                             'parasite proteins reaching them.')
        background = BACKGROUND_FILTERS
        if side == HOST:
            background = st.radio('Test them against', (BACKGROUND_FILTERS,
                                                        BACKGROUND_TISSUES), horizontal=True,
                                  help='Every host protein that came through the filters, '
                                       'or only those expressed where this parasite is. '
                                       'The narrower background asks whether the targets '
                                       'are special among the proteins the parasite can '
                                       'actually meet, and takes the tissues it infects '
                                       'out of the answer; the wider one is the same for '
                                       'every parasite, so two of them can be compared.')
        enrichment = get_enrichment(df_select[df_select['weight'] >= score], data_dir, side,
                                    background, web_utils.get_config_file())
        if not enrichment.empty:
            fdr = st.radio('False discovery rate', (0.01, 0.05, 0.1), horizontal=True,
                           help='The Benjamini-Hochberg corrected significance a process '
                                'has to reach to be counted as enriched.')
            # the figures read the columns the enrichment is built with; only the grid and
            # the file it hands out are renamed for reading
            enrichment_view = enrichment[enrichment['fdr_bh'] <= fdr]
            enrichment_table = enrichment_view[list(ENRICHMENT_COLUMN_NAMES)].rename(
                columns=ENRICHMENT_COLUMN_NAMES)
            st.caption(f'{len(enrichment_table)} processes pass an FDR of {fdr}. Select rows '
                       'to pick them out of the figures below.')
            gb = GridOptionsBuilder.from_dataframe(enrichment_table)
            gb.configure_pagination(paginationAutoPageSize=False,
                                    paginationPageSize=TABLE_PAGE_SIZE)
            gb.configure_grid_options(domLayout='autoHeight')
            gb.configure_side_bar() #Add a sidebar
            gb.configure_selection('multiple', use_checkbox=True, groupSelectsChildren="Group checkbox select children") #Enable multi-row selection
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
        # against the tissues the parasite infects, nothing passing is a result rather than
        # a setting to loosen: the targets are then no more specialised than the proteins
        # the parasite meets anyway, and it is the background that is worth changing
        narrowed = (side == HOST and background == BACKGROUND_TISSUES)
        st.info(f"No biological process passes an FDR of {fdr}. "
                + ("Against the proteins of the tissues this parasite infects, its targets "
                   "carry no process more often than the rest of them. Test against every "
                   "protein that passed the filters to see what its tissues account for."
                   if narrowed else
                   "Loosen the correction above to see the processes that are enriched "
                   "less strongly."))
    elif enrichment_view is not None:
        enrichment_viz = enrichment_view
        if selected_rows is not None and len(selected_rows) > 0:
            selected_terms = selected_rows[GO_TERM_COLUMN].values.tolist()
            enrichment_viz = enrichment_viz[enrichment_viz['go_term'].isin(selected_terms)]

        st.subheader("Enriched biological processes")
        ranked_tab, volcano_tab = st.tabs(["Ranked processes", "All tested processes"])
        with ranked_tab:
            st.caption(f'The {GO_TOP_N} most significantly over-represented biological '
                       f'processes among the {side.lower()} of the network, positioned by '
                       'odds ratio, sized by the number of them annotated to each process '
                       'and shaded by significance. Select rows in the table above to '
                       'restrict the figure.')
            st.plotly_chart(get_enrichment_dotplot(enrichment_viz), width='stretch')
        with volcano_tab:
            st.caption('All processes tested against the network: effect size on the x axis '
                       'and significance on the y axis, with the processes passing the '
                       'selected FDR highlighted.')
            st.plotly_chart(get_enrichment_volcano(enrichment, fdr, selected_terms),
                            width='stretch')

        with st.container():
            if len(selected_terms) > 0:
                if enrichment is not None:
                    highlighted_nodes = enrichment[enrichment['go_term'].isin(selected_terms)]['nodes'].values
                    highlighted_nodes = utils.merge_list_of_lists([i.split(',') for i in highlighted_nodes])
                    highlight_color = {i: HIGHLIGHT_COLOR for i in highlighted_nodes}
                    G = generate_graph(df_select, score,
                                       web_utils.load_protein_annotations(data_dir),
                                       surface_calls)
                    nx.set_node_attributes(G, MUTED_COLOR, 'color')
                    nx.set_node_attributes(G, highlight_color, 'color')
                    # Initiate PyVis network object
                    net = Network(height="450px", width="100%",
                                  bgcolor=NETWORK_BACKGROUND, font_color=LABEL_FONT_COLOR)
                    # Take Networkx graph and translate it to a PyVis graph format
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
        
        fig = get_enrichment_summary(enrichment_view, load_ontology_parents(data_dir))
        st.subheader("Hierarchy of enriched biological processes")
        st.caption('Enriched processes arranged by the Gene Ontology hierarchy: each process '
                   'is nested within the closest enriched process above it, its area is the '
                   'number of network proteins annotated to it and its shade is its '
                   'significance. Click a block to zoom in.')
        st.plotly_chart(fig, width='stretch')

st.markdown("---")
st.markdown("---")

# Footer
with st.container():
    web_utils.footer()