import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils
import web_utils
import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from css import style

st.set_page_config(layout="wide", page_title="OrthoHPI 2.0", menu_items={})
st.session_state.data_dir = 'data'
st.session_state.config_file = 'config.yml'
style.load_css()

web_utils.show_header('Home')

# Read dataset
config = utils.read_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()

# fallback for a parasite without a `group` in the config
UNKNOWN_GROUP = 'Unclassified'
# the parasite groups the surface split can be drawn for. A multicellular parasite reaches
# its host with secreted proteins alone -- the secretome filter keeps it nothing else -- so
# its split is the filter and not the parasite; only the unicellular groups have both
# classes open to them, as every host protein does
UNICELLULAR_GROUPS = ('Apicomplexa', 'Kinetoplastida', 'Other protozoa', 'Microsporidia')
UNKNOWN_COLOR = '#999999'
# A host with two parasites still has to fit two labels under its column, so the columns
# are sized as if every host had at least this many parasites. Human has 35 against the
# two of pig, and strictly proportional columns leave pig a strip its labels overrun.
MIN_COLUMN = 4
# the confidence the counts of the page are drawn at, and the range the slider spans, the
# same default and range as the network page: a parasite counted here then agrees with the
# network the reader opens next instead of being several times larger than it
MIN_SCORE, MAX_SCORE, DEFAULT_SCORE = 0.4, 0.9, 0.4
# share of the height of a figure taken by the strip of taxonomic group under its columns.
# Enough to read as a band of colour, not enough to be read as a quantity of its own
BAND_HEIGHT = 0.06
# and the share of it left blank between the columns and that strip
BAND_GAP = 0.04
# the two surface classes of the secretion figure, as two shades of the blue the page is
# headed in: the bar is one whole split in two, which two shades of a hue say and two
# hues do not. They are also outside the palette of the taxonomic groups in the strip
# under the columns, so neither is read as a clade
SURFACE_COLORS = {'Extracellular': '#a6bddb', 'Cell membrane': '#045a8d'}
# and what to outline a box of that class in, where the class is drawn as a box rather than
# as a bar: the pale shade is a fill and an outline drawn in it on a white background is an
# outline the reader has to look for
SURFACE_LINE_COLORS = {'Extracellular': '#3690c0', 'Cell membrane': '#045a8d'}
# the probability each class was called at, for the dotted line of its column
CLASS_CUTOFFS = {'Extracellular': web_utils.DEEPLOC_EXTRACELLULAR_CUTOFF,
                 'Cell membrane': web_utils.DEEPLOC_MEMBRANE_CUTOFF}


def short_name(parasite):
    '''`Plasmodium vivax` as `P. vivax`, the same abbreviation the circos labels its
    arcs with, so the same parasite reads the same way on every page.'''
    return f'{parasite[0]}. {parasite.split(" ")[1]}'


@st.cache_data(show_spinner=False)
def get_overview_predictions(data_dir, config):
    '''
    Every predicted interaction, labelled with the host species it was predicted against
    and with the taxonomic group of its parasite. A parasite infecting more than one host
    appears once per host, since the point of this page is comparing those.

    :param str data_dir: directory holding predictions.parquet
    :param dict config: parsed configuration
    :return: one row per predicted interaction, with host, parasite, group and weight
    '''
    predictions = web_utils.load_predictions(data_dir)
    groups = {p['label']: p.get('group', UNKNOWN_GROUP) for p in config['parasites'].values()}
    order = {g: i for i, g in enumerate(config.get('parasite_groups', {}))}

    frames = []
    for taxid, host in config['hosts'].items():
        frame = predictions.loc[predictions['taxid2'] == str(taxid),
                                ['taxid1_label', 'weight']]
        if not frame.empty:
            frames.append(frame.assign(host=host['label']))
    df = pd.concat(frames, ignore_index=True)

    df['group'] = df['taxid1_label'].map(lambda p: groups.get(p, UNKNOWN_GROUP))
    # what puts the parasites of a host in the order of the circos: taxonomic group as
    # the configuration declares it, then name, with an unclassified parasite last
    df['group_rank'] = df['group'].map(lambda g: order.get(g, len(order)))
    df['name'] = df['taxid1_label'].map(short_name)

    return df


def host_columns(df, band=False):
    '''
    The skeleton the figures are drawn on: one column per host, as wide as the number of
    parasites infecting it, sharing a y axis so a bar or a box can be compared straight
    across the hosts rather than only within one.

    :param df: overview predictions, as get_overview_predictions builds them
    :param bool band: leave a shallow row under the columns for a strip of the taxonomic
                      group of each parasite, for a figure whose bars carry a colour of
                      their own and cannot also carry the group
    :return: (figure, [(host, its predictions, its parasites in axis order), ...])
    '''
    hosts = []
    for host in df['host'].unique():
        host_df = df[df['host'] == host]
        parasites = host_df[['taxid1_label', 'name', 'group_rank']].drop_duplicates()
        parasites = parasites.sort_values(by=['group_rank', 'taxid1_label'], kind='stable')
        hosts.append((host, host_df, parasites['name'].tolist()))

    widths = [max(len(names), MIN_COLUMN) for _, _, names in hosts]
    figure = make_subplots(rows=2 if band else 1, cols=len(hosts), shared_yaxes=True,
                           column_widths=[w / sum(widths) for w in widths],
                           # the titles fill the top row, which is the row of the figure
                           subplot_titles=[host for host, _, _ in hosts],
                           row_heights=[1 - BAND_HEIGHT, BAND_HEIGHT] if band else None,
                           # enough of a gap that the strip is read as a second thing about
                           # the columns rather than as the foot of the columns themselves
                           vertical_spacing=BAND_GAP, horizontal_spacing=0.02)

    return figure, hosts


def add_group_band(figure, hosts, palette, labelled):
    '''
    The strip of taxonomic group under each column, one segment per parasite in the colour
    the other figures of the page draw that group in. The segments touch, so the strip
    reads as a strip and the clades under a column are blocks rather than a row of bars.

    The groups get a legend of their own (`legend2`), under the one naming the colours of
    the bars: the two say different things about different parts of the figure, and in a
    single row of entries the clades read as more of what the bars are split into.

    :param figure: a figure host_columns built with band=True, modified in place
    :param list hosts: the (host, its rows, its parasites in axis order) of host_columns
    :param dict palette: {taxonomic group: colour}
    :param set labelled: groups already in the legend, added to as they are drawn
    '''
    for column, (host, host_df, names) in enumerate(hosts, start=1):
        for group, rows in host_df.groupby('group', observed=True):
            figure.add_trace(
                go.Bar(x=rows['name'], y=[1] * len(rows), name=group, width=1,
                       marker_color=palette.get(group, UNKNOWN_COLOR),
                       legend='legend2', legendgroup=group,
                       showlegend=group not in labelled,
                       hovertemplate='%{x}' f'<extra>{group}</extra>'),
                row=2, col=column)
            labelled.add(group)
        figure.update_xaxes(categoryorder='array', categoryarray=names, row=2, col=column)
    figure.update_yaxes(visible=False, range=[0, 1], row=2)


def style_host_columns(figure, y_title):
    '''The layout the two figures share: the parasite names under each column, the
    quantity named once down the left, and the taxonomic groups as the legend.'''
    figure.update_layout(height=470, plot_bgcolor='white',
                         margin=dict(l=0, r=0, t=95, b=10),
                         legend=dict(orientation='h', yanchor='bottom', y=1.14, x=0,
                                     title_text='', font=dict(size=11)))
    figure.update_xaxes(tickangle=-60, showgrid=False, tickfont=dict(size=11))
    figure.update_yaxes(showgrid=True, gridcolor='#f0f0f0', zerolinecolor='#e0e0e0')
    figure.update_yaxes(title_text=y_title, row=1, col=1)

    return figure


@st.cache_data(show_spinner=False)
def get_interactor_proteins(data_dir, config, side, score=None):
    '''
    One side of the predicted interactions, protein by protein, with what DeepLoc says
    about where each protein sits: the probability that it is extracellular -- in the space
    the two meet in -- the probability that it is on a cell membrane, and which of the two
    classes those put it in.

    `side` is 'source' for the parasite proteins each parasite reaches its host with, and
    'target' for the host proteins they reach. Both sides went through a localisation
    filter to get here, but not the same one: the parasite side through the secretome
    filter, which allows a multicellular parasite nothing but secreted proteins, and the
    host side through apply_deeploc_filter, which allows every host both classes. Only the
    host side can therefore be read as a comparison between parasites.

    One row per protein and not per interaction: a parasite protein reaching eleven host
    proteins is one protein, not eleven. A parasite infecting two hosts has its proteins
    counted under each, since the page compares the hosts.

    :param str data_dir: directory holding predictions.parquet and the localisations
    :param dict config: parsed configuration
    :param str side: 'source' for the parasite proteins, 'target' for the host proteins
    :param score: keep only the proteins of interactions predicted at or above this
                  confidence, or None for every prediction. The figures of the page take
                  it both ways: the proportions are read at the threshold, the boxes of
                  the probability itself at every prediction, since a box standing on the
                  eight proteins a threshold leaves is not a distribution
    :return: one row per host, parasite and protein, or None without a localisation table
    '''
    localisations = web_utils.load_deeploc_localisations(data_dir)
    if localisations.empty:
        return None

    predictions = web_utils.load_predictions(data_dir)
    if score is not None:
        predictions = predictions[predictions['weight'] >= score]
    surface = localisations.assign(surface=web_utils.classify_surface(localisations))

    frames = []
    for taxid, host in config['hosts'].items():
        frame = predictions.loc[predictions['taxid2'] == str(taxid),
                                ['taxid1_label', side]].drop_duplicates()
        if not frame.empty:
            frames.append(frame.assign(host=host['label']))
    df = pd.concat(frames, ignore_index=True).merge(surface, left_on=side,
                                                    right_on='protein', how='inner')

    groups = {p['label']: p.get('group', UNKNOWN_GROUP) for p in config['parasites'].values()}
    order = {g: i for i, g in enumerate(config.get('parasite_groups', {}))}
    df['group'] = df['taxid1_label'].map(lambda p: groups.get(p, UNKNOWN_GROUP))
    df['group_rank'] = df['group'].map(lambda g: order.get(g, len(order)))
    df['name'] = df['taxid1_label'].map(short_name)

    return df


@st.cache_data(show_spinner=False)
def get_surface_counts(proteins, every=None):
    '''
    How many of each parasite's proteins fall in each of the surface classes.

    The figure splits its columns over the two surface classes alone, so a parasite is
    measured on what was called rather than on how much of the proteome the model was sure
    about. The proteins in neither class are counted here all the same and are read in the
    hover; a parasite with nothing in either has no column, which on the host side happens
    to nobody -- every host protein is here because it was called one or the other.

    :param proteins: the proteins of one side, as get_interactor_proteins builds them
    :param every: the same proteins before the confidence threshold, to keep a parasite the
                  threshold emptied as a column with nothing in it. Dropping it instead
                  would take a column out of this figure and leave it in the figures above
                  and below, which are read across the page as the same columns
    :return: one row per host and parasite
    '''
    counts = proteins.pivot_table(index=['host', 'taxid1_label', 'name', 'group',
                                         'group_rank'],
                                  columns='surface', values='protein', aggfunc='count',
                                  fill_value=0)
    for surface_class in [web_utils.EXTRACELLULAR, web_utils.CELL_MEMBRANE,
                          web_utils.NOT_SURFACE]:
        if surface_class not in counts.columns:
            counts[surface_class] = 0
    counts = counts.reset_index()

    if every is not None:
        keys = ['host', 'taxid1_label', 'name', 'group', 'group_rank']
        counts = (every[keys].drop_duplicates()
                  .merge(counts, on=keys, how='left').fillna(0))

    return counts


@st.cache_data(show_spinner=False)
def generate_surface_split_per_parasite(df, palette, y_title='host proteins reached',
                                       hover_noun='the host proteins it reaches'):
    '''
    How a parasite's proteins are split between the two surface classes, in the same columns
    and the same order as the figures around it, so they are read together. Drawn for either
    side: the host proteins a parasite reaches, or the proteins of the parasite itself.

    :param df: surface counts, as get_surface_counts builds them
    :param dict palette: {taxonomic group: colour} for the strip under the columns
    :param str y_title: what the columns are a proportion of, named down the left
    :param str hover_noun: the same, phrased for the hover of a bar

    Every column is the whole of what that parasite reaches on its host and is split in two
    by where DeepLoc puts those proteins, so the columns are compared on the split itself
    rather than on how many proteins a parasite reaches: a column that is nearly solid dark
    is a parasite that docks onto the surface of the host cell, and the pale part of one is
    what it meets in the matrix and the fluid around the cell instead. Both classes were
    open to every host protein -- apply_deeploc_filter keeps a host protein for either --
    so the difference between the columns is a difference between the parasites.

    The colour is spent on the two classes, so the taxonomic group each parasite belongs to
    moves to the strip under the columns and the clades of a host are read there as blocks
    of colour.
    '''
    figure, hosts = host_columns(df, band=True)
    labelled = set()
    for column, (host, host_df, names) in enumerate(hosts, start=1):
        counted = host_df[web_utils.EXTRACELLULAR] + host_df[web_utils.CELL_MEMBRANE]
        for surface_class in SURFACE_COLORS:
            figure.add_trace(
                go.Bar(x=host_df['name'], y=host_df[surface_class] / counted,
                       name=surface_class, marker_color=SURFACE_COLORS[surface_class],
                       customdata=host_df[[surface_class]],
                       legendgroup=surface_class, showlegend=surface_class not in labelled,
                       hovertemplate='%{x}<br>%{y:.0%} of ' + hover_noun +
                                     ' (%{customdata[0]} of them)'
                                     f'<extra>{surface_class}</extra>'),
                row=1, col=column)
            labelled.add(surface_class)
        figure.update_xaxes(categoryorder='array', categoryarray=names, row=1, col=column)
    add_group_band(figure, hosts, palette, labelled)

    figure = style_host_columns(figure, y_title)
    # two legends, one above the other and each named, since the colours of the bars and
    # the colours of the strip are two different keys to two different parts of the figure
    figure.update_layout(barmode='stack', bargap=0.2, margin=dict(t=125),
                         legend=dict(y=1.28, title_text='DeepLoc',
                                     title_font=dict(size=11)),
                         legend2=dict(orientation='h', yanchor='bottom', y=1.14, x=0,
                                      title_text='taxonomic group',
                                      title_font=dict(size=11), font=dict(size=11)))
    figure.update_yaxes(range=[0, 1], tickformat='.0%', row=1)
    # the names belong under the strip, which is the foot of the figure now, and a column
    # labelled twice is a column labelled once too often. The room they need is taken off
    # the figure rather than out of the margin, which is where the strip has pushed them
    figure.update_xaxes(showticklabels=False, row=1)
    figure.update_xaxes(automargin=True, row=2)
    # and the same for the percentages down the left, which the two figures above have no
    # room for either -- there they are a count anyone can read off the bars, here they
    # are the scale the split is read on
    figure.update_yaxes(automargin=True, row=1, col=1)

    return figure


@st.cache_data(show_spinner=False)
def generate_surface_scores_per_parasite(proteins, palette, score, cutoff, y_title,
                                         point_size=2.5):
    '''
    The spread of the probability itself, before it is a class: one box per parasite over
    the proteins of one surface class, in the same columns and colours as the figures above
    it. Read within a class rather than across the whole scale, a box says how sure DeepLoc
    was of the proteins it did call that -- a box sitting on the cut-off is a parasite whose
    proteins only just qualified, a box near 1 one whose proteins are unambiguous.

    A parasite with no protein of the class has no box, which is how the membrane figure
    comes out as the unicellular parasites alone: the secretome filter keeps a multicellular
    parasite's proteins for being secreted and for nothing else.

    Every protein is drawn as well as summarised, jittered across the width of its box.
    Several of these boxes stand on ten or twenty proteins, which is too few for a box to
    be read as a distribution without seeing what is behind it. The confidence figure at
    the top of the page keeps its points off for the opposite reason: there a point is an
    interaction rather than a protein, and there are thousands of them to a column.

    The scale starts a little below the cut-off, since nothing of the class is under it and
    a box drawn against the whole 0 to 1 is a box flattened into the top of the figure.

    :param proteins: the proteins of one surface class, as get_parasite_proteins builds them
    :param dict palette: {taxonomic group: colour}
    :param str score: the probability column to draw
    :param float cutoff: the cut-off of that class, drawn as a line and the foot of the scale
    :param str y_title: what the probability is called down the left of the figure
    :param float point_size: diameter of a protein drawn beside its box, in pixels. The
                             proteins are jittered across the width of the box and the box
                             itself is left unfilled so they are read through it; a column
                             narrower than its proteins are many comes out as a cloud rather
                             than as points that can be counted
    '''
    figure, hosts = host_columns(proteins)
    labelled = set()
    for column, (host, host_df, names) in enumerate(hosts, start=1):
        for group, rows in host_df.groupby('group', observed=True):
            figure.add_trace(
                go.Box(x=rows['name'], y=rows[score], name=group,
                       marker=dict(color=palette.get(group, UNKNOWN_COLOR),
                                   size=point_size, opacity=0.45),
                       line=dict(color=palette.get(group, UNKNOWN_COLOR), width=1),
                       fillcolor='rgba(0,0,0,0)',
                       boxpoints='all', jitter=0.8, pointpos=0,
                       legendgroup=group, showlegend=group not in labelled,
                       hovertemplate='%{x}<br>' f'{y_title} ' '%{y:.2f}'
                                     f'<extra>{host}</extra>'),
                row=1, col=column)
            labelled.add(group)
        figure.update_xaxes(categoryorder='array', categoryarray=names, row=1, col=column)

    figure = style_host_columns(figure, y_title)
    figure.update_yaxes(range=[cutoff - 0.05, 1], automargin=True, row=1, col=1)
    # the names are as long as the columns are narrow, so the room they need comes off the
    # figure rather than out of the margin they would otherwise be cut in
    figure.update_xaxes(automargin=True)
    figure.add_hline(y=cutoff, line_width=1, line_dash='dot', line_color='#969696')

    return figure


@st.cache_data(show_spinner=False)
def generate_host_score_boxes(proteins, point_size=3):
    '''
    How sure DeepLoc was of the host proteins, one column per surface class and one box per
    host inside it.

    Grouped by class and not by host: the two classes are called on two different
    probabilities and at two different cut-offs, so what a box means changes with the
    class and not with the host. A column per class puts the hosts on a common axis and
    leaves one cut-off to draw per column instead of both in every column.

    Per host and not per parasite. The parasites of a host draw their interactors from the
    same few hundred host proteins, so a box per parasite is a box over a sample of one
    pool: their medians span 0.65 to 0.79 with no order to them, forty-six near-copies of
    the pool and of each other.

    A host protein is counted once however many parasites reach it, since here it is a
    protein of the host and not an interactor of anything.

    :param proteins: the host proteins, as get_interactor_proteins(side='target') builds them
    :param float point_size: diameter of a protein drawn beside its box, in pixels
    '''
    scored = proteins.drop_duplicates(['host', 'protein']).copy()
    scored = scored[scored['surface'] != web_utils.NOT_SURFACE]
    scored['score'] = np.where(scored['surface'] == web_utils.EXTRACELLULAR,
                               scored['extracellular'], scored['cell_membrane'])
    # host_columns splits on 'host' and orders what is in a column by 'group_rank'; here a
    # column is a surface class and what is in it are the hosts, so the two are swapped
    host_order = {host: rank for rank, host in enumerate(scored['host'].unique())}
    scored['name'] = scored['host']
    scored['taxid1_label'] = scored['host']
    scored['group_rank'] = scored['host'].map(host_order)
    scored['host'] = scored['surface']

    figure, classes = host_columns(scored)
    for column, (surface_class, class_df, names) in enumerate(classes, start=1):
        figure.add_trace(
            go.Box(x=class_df['name'], y=class_df['score'], name=surface_class,
                   marker=dict(color=SURFACE_COLORS[surface_class], size=point_size,
                               opacity=0.45),
                   line=dict(color=SURFACE_LINE_COLORS[surface_class], width=1),
                   fillcolor='rgba(0,0,0,0)', boxpoints='all', jitter=0.8, pointpos=0,
                   showlegend=False,
                   hovertemplate='%{x}<br>P %{y:.2f}' f'<extra>{surface_class}</extra>'),
            row=1, col=column)
        figure.update_xaxes(categoryorder='array', categoryarray=names, row=1, col=column)
        # the cut-off this class was called at, drawn in the column it applies to rather
        # than across the figure, where it would read as a threshold on both
        figure.add_hline(y=CLASS_CUTOFFS[surface_class], line_width=1, line_dash='dot',
                         line_color=SURFACE_LINE_COLORS[surface_class], row=1, col=column)

    figure = style_host_columns(figure, 'P(the class it was called)')
    figure.update_yaxes(range=[min(CLASS_CUTOFFS.values()) - 0.05, 1], automargin=True,
                        row=1, col=1)
    figure.update_xaxes(tickangle=0, automargin=True)

    return figure


@st.cache_data(show_spinner=False)
def generate_interactions_per_parasite(df, palette, score):
    '''
    How many interactions are predicted for each parasite at or above a confidence, in one
    column per host. Bars are coloured by taxonomic group.

    The columns are laid out from every prediction and only the bars are thresholded, so a
    parasite left with nothing keeps its place on the axis and is read as an absence rather
    than disappearing out of a figure the one below it still draws it in.

    :param df: overview predictions, as get_overview_predictions builds them
    :param dict palette: {taxonomic group: colour}
    :param float score: the confidence to count from
    '''
    figure, hosts = host_columns(df)
    labelled_groups = set()
    for column, (host, host_df, names) in enumerate(hosts, start=1):
        kept = host_df[host_df['weight'] >= score]
        counts = kept.groupby(['name', 'group'], observed=True).size().reset_index(name='count')
        for group, rows in counts.groupby('group', observed=True):
            figure.add_trace(
                go.Bar(x=rows['name'], y=rows['count'], name=group,
                       marker_color=palette.get(group, UNKNOWN_COLOR),
                       legendgroup=group, showlegend=group not in labelled_groups,
                       hovertemplate='%{x}<br>%{y} predicted interactions'
                                     f'<extra>{group}</extra>'),
                row=1, col=column)
            labelled_groups.add(group)
        figure.update_xaxes(categoryorder='array', categoryarray=names, row=1, col=column)

    return style_host_columns(figure, 'predicted interactions')


@st.cache_data(show_spinner=False)
def generate_confidence_per_parasite(df, palette, score):
    '''
    The spread of the confidence score of each parasite's predicted interactions, in the
    same columns and colours as the counts above, so the two figures are read together:
    a parasite with many interactions and a low box has many weakly supported ones.

    The one figure of the page the slider does not filter, drawn on every prediction with
    the slider as a line across it instead. Filtering a distribution by the axis it is
    drawn on says nothing -- every lower whisker becomes the slider -- where the line makes
    this figure the key to the slider: how much of a parasite's box stands above it is how
    much of that parasite survives into the counts above.

    The outlying points are left off. Every interaction is a point, and a few hundred of
    them beside each box hide the boxes themselves.

    :param df: overview predictions, as get_overview_predictions builds them
    :param dict palette: {taxonomic group: colour}
    :param float score: where the counts above are cut, drawn as a line
    '''
    figure, hosts = host_columns(df)
    labelled = set()
    for column, (host, host_df, names) in enumerate(hosts, start=1):
        for group, rows in host_df.groupby('group', observed=True):
            figure.add_trace(
                go.Box(x=rows['name'], y=rows['weight'], name=group,
                       marker_color=palette.get(group, UNKNOWN_COLOR),
                       line_width=1, boxpoints=False,
                       legendgroup=group, showlegend=group not in labelled,
                       hovertemplate='%{x}<br>score %{y}' f'<extra>{host}</extra>'),
                row=1, col=column)
            labelled.add(group)
        figure.update_xaxes(categoryorder='array', categoryarray=names, row=1, col=column)

    figure = style_host_columns(figure, 'confidence score')
    figure.add_hline(y=score, line_width=1, line_dash='dot', line_color='#969696',
                     annotation_text='counted above', annotation_position='top left',
                     annotation_font=dict(size=10, color='#969696'))

    return figure


st.caption('Protein-protein interactions between parasites and their hosts, predicted by '
           'orthology transfer and restricted to the host proteins expressed in a tissue '
           'the parasite is known to infect. This page presents every prediction per host; '
           '**Parasites of a host** compares the parasites of a single host, **Hosts of a '
           'parasite** follows one parasite across the hosts it infects, and **Host-parasite '
           'network** shows the network of one host-parasite pair.')
st.markdown("---")

overview = get_overview_predictions(data_dir, config)
parasite_palette = config.get('parasite_groups', {})

# one slider for the counts and the proportions of the page, so what is counted here is
# what the network page opens on rather than several times more of it. It sits in the
# middle of three columns, as the sliders of the other pages do: left to itself a slider
# takes the whole width of the page, which is a metre of track for a range of half a point
with st.columns(3)[1]:
    score = st.slider('Confidence score', MIN_SCORE, MAX_SCORE, DEFAULT_SCORE,
                      help='Interactions predicted below this confidence are left out of '
                           'the counts and the proportions. The boxplots keep every '
                           'prediction: the confidence figure draws this threshold as a '
                           'line instead, and the localization figures stand on too few '
                           'proteins to be thresholded as well.')

st.subheader("Number of predicted interactions per parasite")
st.caption('Predicted interactions per parasite at or above the confidence set above, grouped '
           'by host and coloured by parasite taxonomic group.')
st.plotly_chart(generate_interactions_per_parasite(overview, parasite_palette, score),
                width='stretch')

st.subheader("Confidence of the predicted interactions per parasite")
st.caption('Boxplots of the distribution of confidence scores per parasite. Scores derive from '
           'the evidence supporting the orthologous interaction from which each prediction was '
           'transferred. Every prediction is counted here, whatever the slider is set to; the '
           'dotted line marks it, so the part of a box above the line is the part of that '
           'parasite counted in the figure above.')
st.plotly_chart(generate_confidence_per_parasite(overview, parasite_palette, score),
                width='stretch')

# the localisation figures are drawn twice over: the proportions at the threshold, and the
# boxes of the probability itself on every prediction. At 0.7 a parasite is left a median
# of seventeen host proteins and eight of its own, which a proportion can still be read
# off and a box cannot
host_proteins = get_interactor_proteins(data_dir, config, 'target')
parasite_proteins = get_interactor_proteins(data_dir, config, 'source')
host_proteins_kept = get_interactor_proteins(data_dir, config, 'target', score)
parasite_proteins_kept = get_interactor_proteins(data_dir, config, 'source', score)
if host_proteins is not None:
    st.subheader("Proportion of extracellular and membrane host proteins")
    st.caption('Subcellular localization predicted by DeepLoc 2 for the host proteins each '
               'parasite reaches at or above the confidence set above, divided into cell '
               'membrane and extracellular. The strip below the columns indicates taxonomic '
               'group, coloured as above.')
    st.plotly_chart(
        generate_surface_split_per_parasite(get_surface_counts(host_proteins_kept,
                                                               every=host_proteins),
                                            parasite_palette), width='stretch')

    st.subheader("Localization confidence of host proteins")
    st.caption('Boxplots of the DeepLoc 2 probabilities of the host proteins for the class '
               'they were assigned to, one column per class and one box per host. The dotted '
               'line in each column marks the cut-off that class was called at. Each host '
               'protein is counted once, irrespective of the number of parasites reaching it, '
               'and every prediction is counted whatever the slider is set to.')
    st.plotly_chart(generate_host_score_boxes(host_proteins), width='stretch')

if parasite_proteins is not None:
    unicellular = parasite_proteins[parasite_proteins['group'].isin(UNICELLULAR_GROUPS)]
    kept_unicellular = parasite_proteins_kept[
        parasite_proteins_kept['group'].isin(UNICELLULAR_GROUPS)]
    if not unicellular.empty:
        st.subheader("Proportion of extracellular and membrane parasite proteins")
        st.caption('Subcellular localization predicted by DeepLoc 2 for the proteins each '
                   'unicellular parasite reaches its host with at or above the confidence set '
                   'above, divided into cell membrane and extracellular. Multicellular '
                   'parasites are omitted, as the secretome filter admits only their secreted '
                   'proteins. The strip below the columns indicates taxonomic group, coloured '
                   'as above.')
        st.plotly_chart(
            generate_surface_split_per_parasite(get_surface_counts(kept_unicellular,
                                                                   every=unicellular),
                                                parasite_palette,
                                                y_title='proteins of the parasite',
                                                hover_noun='its proteins'),
            width='stretch')

    st.subheader("Localization confidence of extracellular parasite proteins")
    st.caption('Boxplots of the DeepLoc 2 probability of extracellular localization for the '
               'extracellular proteins of each parasite, over every prediction whatever the '
               'slider is set to. The dotted line marks the cut-off. Individual proteins are '
               'shown as points behind each box; for parasites with a hundred or more proteins '
               'the points are read as density.')
    st.plotly_chart(
        generate_surface_scores_per_parasite(
            parasite_proteins[parasite_proteins['surface'] == web_utils.EXTRACELLULAR],
            parasite_palette, 'extracellular', web_utils.DEEPLOC_EXTRACELLULAR_CUTOFF,
            # forty-five columns to a row and up to a hundred and eighty proteins in one of
            # them, so the smallest dot that still carries colour
            'P(extracellular)', point_size=2),
        width='stretch')

    st.subheader("Localization confidence of membrane parasite proteins")
    st.caption('The equivalent for the parasite proteins assigned to cell membrane, scored on '
               'that probability, again over every prediction. Only unicellular parasites are '
               'represented. Individual proteins are shown as points; boxes over very few proteins (six for *C. '
               'parvum*, one for *V. corneae*) should not be read as distributions.')
    st.plotly_chart(
        generate_surface_scores_per_parasite(
            parasite_proteins[parasite_proteins['surface'] == web_utils.CELL_MEMBRANE],
            parasite_palette, 'cell_membrane', web_utils.DEEPLOC_MEMBRANE_CUTOFF,
            # a third of the parasites and a fifth of the proteins of the figure above, so
            # the columns are wide enough for the proteins to be told apart
            'P(cell membrane)', point_size=3),
        width='stretch')

st.markdown("---")


# Footer
with st.container():
    web_utils.footer()
