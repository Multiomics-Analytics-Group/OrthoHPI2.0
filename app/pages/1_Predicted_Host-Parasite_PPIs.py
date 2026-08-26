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
# row. The host species is put back when the selected host covers more than one (Rodent is
# rat and mouse), where it is what tells the rows apart.
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
    :param str host: the selected host group
    :return: the table with the two columns in front, or unchanged if it is empty
    '''
    if table.empty:
        return table

    # a host group covering two species already carries the species of each row, and that
    # is the column the group name must not overwrite: it is the one thing the rows of a
    # pooled host differ by
    host_column = 'Host group' if 'Host species' in table.columns else 'Host species'
    named = table.assign(**{'Parasite species': parasite, host_column: host})
    front = ['Parasite species', host_column]

    return named[front + [c for c in named.columns if c not in front]]


def generate_tissue_filters(df):
    options = df['Tissue'].unique().tolist()

    return options

def generate_cell_type_filters(df):
    options = df['Cell type'].dropna().unique().tolist()

    return options

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

@st.cache_data(max_entries=3, ttl=1800)
def get_enrichment(pred_df, data_dir):
    species = pred_df['taxid1'].unique().tolist() + pred_df['taxid2'].unique().tolist()
    species = [int(s) for s in species]
    # the filter is pushed down to the reader (fastparquet prunes row groups only,
    # so the exact selection is still applied afterwards)
    go_df = utils.read_parquet_file(input_file=f'{data_dir}/gos.parquet', filters=[('taxid', 'in', species)])
    go_df = go_df[go_df['taxid'].isin(species)]
    enrichment = utils.calculate_enrichment(pred_df, go_df)
    # A is the number of proteins of the network annotated to the term, which is what
    # every figure below sizes its marks by
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
    '''child term -> its parent terms, the shape the ontology is walked upward in.'''
    ontology = load_ontology(data_dir)

    return ontology.groupby('child')['parent'].apply(list).to_dict()


def nearest_enriched_ancestors(terms, parents_of):
    '''
    For each enriched term, the closest term above it in the ontology that is also
    enriched, so the processes can be nested in each other. Terms with no enriched
    ancestor sit at the top level. The ontology is a graph rather than a tree -- a term
    can have several parents -- so it is walked breadth-first and the first enriched
    level reached wins.
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
    view = prepare_enrichment_view(enrichment_df).drop_duplicates(subset='go_term')
    terms = view['go_term'].tolist()
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
        ids=[GO_TREEMAP_ROOT_ID] + terms,
        # unwrapped: a block draws its name on one line and drops it when it does not
        # fit, which keeps the header of a group readable instead of hiding it
        labels=[GO_TREEMAP_ROOT_LABEL] + terms,
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

    :param list edges: vis.js edge dictionaries, as pyvis built them
    :param df: predictions dataframe of the selected parasite
    :param dict annotations: STRING id --> descriptive protein name
    :return: the same edges, each with the proteins of its interaction added
    '''
    # the species travels with the edge as well, so the dialog can name it beside each
    # protein and say which of the two models is the parasite's and which the host's
    cols = ['source', 'source_name', 'source_uniprot', 'taxid1_label',
            'target', 'target_name', 'target_uniprot', 'taxid2_label', 'weight']
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
            'weight': float(row.weight)}

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
    clicked in the network. Opened over the network rather than placed under it so that
    the table and the enrichment plots below do not move on every click.

    :param dict edge: the clicked edge, as annotated by annotate_edges
    '''
    st.markdown(f"**{edge['parasite']}** ({edge['parasite_species']}) &ndash; "
                f"**{edge['host']}** ({edge['host_species']}) &nbsp;·&nbsp; "
                f"interaction confidence score {edge['weight']:.2f}", unsafe_allow_html=True)
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
        score = st.slider('Confidence score', 0.4, 0.9, 0.7)

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
        cell_type_options = generate_cell_type_filters(df_select)
        if len(cell_type_options) > 0:
            selected_cell_types = st.multiselect(
                'Select cell types to filter the predicted PPI', cell_type_options,
                help='The single cell annotation is human, so this is offered for human '
                     'host proteins alone. A host protein with no cell type annotation is '
                     'left out once a cell type is chosen.')
            if len(selected_cell_types) > 0:
                df_select = df_select[df_select['Cell type'].isin(selected_cell_types)]

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
                                     selected_taxids)


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
        gb.configure_pagination(paginationAutoPageSize=True) #Add pagination
        gb.configure_side_bar() #Add a sidebar
        gridOptions = gb.build()
        grid_response = AgGrid(
                            table,
                            gridOptions=gridOptions,
                            data_return_mode='AS_INPUT',
                            fit_columns_on_grid_load=False,
                            enable_enterprise_modules=True,
                            height=350
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
        st.caption('Biological processes over-represented among the proteins of the network, '
                   'host and parasite alike, against every annotated protein of the two '
                   "species. Fisher's exact test, corrected across terms with "
                   'Benjamini-Hochberg; only processes carried by 11 to 499 proteins of the '
                   'network are tested.')
        enrichment = get_enrichment(df_select[df_select['weight'] >= score], data_dir)
        if not enrichment.empty:
            fdr = st.radio('False discovery rate', (0.01, 0.05, 0.1), horizontal=True,
                           help='The Benjamini-Hochberg corrected significance a process '
                                'has to reach to be counted as enriched.')
            enrichment_table = enrichment[enrichment['fdr_bh'] <= fdr][
                list(ENRICHMENT_COLUMN_NAMES)].rename(columns=ENRICHMENT_COLUMN_NAMES)
            st.caption(f'{len(enrichment_table)} processes pass an FDR of {fdr}. Select rows '
                       'to pick them out of the figures below.')
            gb = GridOptionsBuilder.from_dataframe(enrichment_table)
            gb.configure_pagination(paginationAutoPageSize=True) #Add pagination
            gb.configure_side_bar() #Add a sidebar
            gb.configure_selection('multiple', use_checkbox=True, groupSelectsChildren="Group checkbox select children") #Enable multi-row selection
            gridOptions = gb.build()
            grid_response = AgGrid(
                                enrichment_table,
                                gridOptions=gridOptions,
                                data_return_mode='AS_INPUT',
                                fit_columns_on_grid_load=False,
                                enable_enterprise_modules=True,
                                height=350
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
    if enrichment_table is not None and enrichment_table.empty:
        st.info(f"No biological process passes an FDR of {fdr}. Loosen the correction above "
                "to see the processes that are enriched less strongly.")
    elif enrichment_table is not None:
        enrichment_viz = enrichment_table
        if selected_rows is not None and len(selected_rows) > 0:
            selected_terms = selected_rows[GO_TERM_COLUMN].values.tolist()
            enrichment_viz = enrichment_viz[enrichment_viz['go_term'].isin(selected_terms)]

        st.subheader("Enriched biological processes")
        ranked_tab, volcano_tab = st.tabs(["Ranked processes", "All tested processes"])
        with ranked_tab:
            st.caption(f'The {GO_TOP_N} most significantly over-represented biological '
                       'processes among the host proteins of the network, positioned by odds '
                       'ratio, sized by the number of network proteins annotated to them and '
                       'shaded by significance. Select rows in the table above to restrict '
                       'the figure.')
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
        
        fig = get_enrichment_summary(enrichment_table, load_ontology_parents(data_dir))
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