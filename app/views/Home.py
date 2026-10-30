import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import utils
import web_utils
import hosts_parasites
import similarity_heatmap
import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

web_utils.show_header()

config = utils.read_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()

# fallback for a parasite without a `group` in the config
UNKNOWN_GROUP = 'Unclassified'
# the membrane figure is drawn over these: a multicellular parasite reaches its host with
# secreted proteins alone
UNICELLULAR_GROUPS = ('Apicomplexa', 'Kinetoplastida', 'Metamonada', 'Amoebozoa',
                      'Microsporidia')
UNKNOWN_COLOR = '#999999'
# columns are sized as if every host had at least this many parasites, so a two-parasite
# host still fits its labels
MIN_COLUMN = 4
# caption of the niche strips under the columns
INTRACELLULAR_NOTE = ('A parasite with an intracellular stage in the host counts as '
                      'intracellular, since the reach of that stage is the wider of the two.')
NICHE_STRIP = ('The strip below the columns indicates whether the parasite lives inside a '
               'host cell or outside one. ' + INTRACELLULAR_NOTE)
# and of the counts figure, which carries the taxonomic group as well
BANDS_STRIP = ('The upper strip below the columns indicates the taxonomic group of the '
               'parasite, the lower one whether it lives inside a host cell or outside '
               'one. ' + INTRACELLULAR_NOTE)
# subplot titles are shrunk and wrapped until they fit their column; about 0.55 of the font
# size a character
TITLE_CHAR = 0.55
# pixels of a column its title is not written across
COLUMN_PADDING = 12
# plotly's own size for a subplot title
TITLE_SIZE = 16
# below this a host name is no longer read, so the title is drawn wider than the column
# instead
SMALLEST_TITLE = 10
# share of the figure height taken by the strip under its columns
BAND_HEIGHT = 0.06
# and the share left blank between the columns and that strip
BAND_GAP = 0.04
SURFACE_COLORS = web_utils.LOCALISATION_COLORS
# outline of a box of that class, the pale fills being invisible as outlines; no entry for
# the mixed classes
SURFACE_LINE_COLORS = {'Extracellular': '#3690c0', 'Cell membrane': '#045a8d',
                       'Cytoplasm': '#e6550d', 'Nucleus': '#a63603'}
# the host side carries the four classes its filter reads, the parasite side the surface
# pair the secretome filter selected on
HOST_SPLIT_CLASSES = web_utils.HOST_CLASSES + (web_utils.SEVERAL,)
PARASITE_SPLIT_CLASSES = web_utils.SURFACE_CLASSES + (web_utils.BOTH_SURFACE,)
# what the counts figure is coloured by: the taxonomic group, or the localization of one end
# of each interaction
SPLIT_SIDES = {'Taxonomic group': (None, None),
               'Host proteins': ('target', HOST_SPLIT_CLASSES),
               'Parasite proteins': ('source', PARASITE_SPLIT_CLASSES)}
# how far under the lowest cut-off or point a probability scale starts
SCALE_MARGIN = 0.05
# rows of the shared orthologous group dot plot, pixels a row takes, and what the strip,
# the names and the legends take besides
TOP_OGS = 40
DOT_ROW = 19
DOT_CHROME = 230
# characters a line of an orthologous group name is wrapped at
LABEL_CHARS = 30
# room a group name takes left of the dots: ~6.2 px a character, plus the axis title
LABEL_CHAR = 6.2
LABEL_PADDING = 40
# height of one horizontal legend row, in pixels
LEGEND_ROW = 46
# session key counting remounts of the cross-host overlap matrix, so the same cell can be
# opened twice running
CELL_NONCE_KEY = 'cross_host_cell_nonce'
# side of that matrix in pixels at most
OVERLAP_SIDE = 640
# what the rows of that matrix can be ordered on, first to last, and the sort keys of each;
# the parasite and its hosts close every order, so a parasite or genus stays together
# outside the host order
OVERLAP_ORDERS = {
    web_utils.NICHE_TITLE: ['niche_rank', 'group_rank', 'taxid1_label', 'host_rank'],
    'taxonomic group': ['group_rank', 'niche_rank', 'taxid1_label', 'host_rank'],
    'host': ['host_rank', 'niche_rank', 'group_rank', 'taxid1_label']}
# the five numbers a box is drawn from
BOX_QUANTILES = [0, 0.25, 0.5, 0.75, 1]


def score_floor(*values):
    '''
    Where a probability scale starts: a little under the lowest of what has to be
    visible on it.
    '''
    return max(0.0, min(values) - SCALE_MARGIN)


def short_name(parasite):
    '''
    `Plasmodium vivax` as `P. vivax`, the same abbreviation the circos labels its arcs
    with, so the same parasite reads the same way on every page.
    '''
    return f'{parasite[0]}. {parasite.split(" ")[1]}'


def fit_titles(hosts, rooms):
    '''The names of the hosts, written to fit the columns they head.'''
    def forms(host):
        species, _, common = host.partition(' (')
        common = f'({common}' if common else ''
        if not common:
            return [[species]] * 3

        return [[f'{species} {common}'], [species, common], [short_name(species), common]]

    def size_of(lines, room):
        return min(TITLE_SIZE, int(room / (TITLE_CHAR * max(len(line) for line in lines))))

    written = [[forms(host)[form] for host in hosts] for form in range(3)]
    # the size a form can be written at is what its narrowest column allows; the earliest
    # form wins a tie
    sizes = [min(size_of(lines, room) for lines, room in zip(form, rooms))
             for form in written]
    best = sizes.index(max(sizes))

    return ['<br>'.join(lines) for lines in written[best]], max(sizes[best], SMALLEST_TITLE)


@st.cache_data(show_spinner=False)
def get_overview_predictions(data_dir, config):
    '''
    Every predicted interaction, labelled with the host species it was predicted against
    and with the taxonomic group of its parasite.
    '''
    predictions = web_utils.load_predictions(data_dir)
    groups = {p['label']: p.get('group', UNKNOWN_GROUP) for p in config['parasites'].values()}
    order = {g: i for i, g in enumerate(config.get('parasite_groups', {}))}

    frames = []
    for taxid, host in config['hosts'].items():
        frame = predictions.loc[predictions['taxid2'] == str(taxid),
                                ['taxid1_label', 'weight', 'source', 'target', 'target_name',
                                 'group2']]
        if not frame.empty:
            frames.append(frame.assign(host=host['label']))
    df = pd.concat(frames, ignore_index=True)

    df['group'] = df['taxid1_label'].map(lambda p: groups.get(p, UNKNOWN_GROUP))
    # the order of the paper's figures: niche, then taxonomic group as configured, then
    # name, unclassified last
    df['group_rank'] = df['group'].map(lambda g: order.get(g, len(order)))
    df['name'] = df['taxid1_label'].map(short_name)
    df['niche'] = df['taxid1_label'].map(web_utils.get_niches(config)).fillna(
        web_utils.UNKNOWN_NICHE)

    return df


def host_columns(df, width, bands=0, band_height=BAND_HEIGHT):
    '''
    The skeleton the figures are drawn on: one column per host, as wide as the number of
    parasites infecting it, sharing a y axis so a bar or a box can be compared straight
    across the hosts rather than only within one. `band_height` is the share of the
    figure a strip under the columns takes, which a tall figure hands over in pixels.
    '''
    hosts = []
    for host in df['host'].unique():
        host_df = df[df['host'] == host]
        parasites = host_df[['taxid1_label', 'name', 'group_rank']].drop_duplicates(
            'taxid1_label')
        parasites['niche_rank'] = parasites.index.map(host_df['niche']).map(web_utils.niche_rank)
        parasites = parasites.sort_values(by=['niche_rank', 'group_rank', 'taxid1_label'],
                                          kind='stable')
        hosts.append((host, host_df, parasites['name'].tolist()))

    widths = [max(len(names), MIN_COLUMN) for _, _, names in hosts]
    figure = make_subplots(rows=1 + bands, cols=len(hosts), shared_yaxes=True,
                           column_widths=[w / sum(widths) for w in widths],
                           subplot_titles=[host for host, _, _ in hosts],
                           row_heights=([1 - bands * band_height] + [band_height] * bands
                                        if bands else None),
                           # enough of a gap that a strip is read as a second thing about
                           # the columns
                           vertical_spacing=BAND_GAP, horizontal_spacing=0.02)

    # the subplot titles are the first annotations, one per column; refitted here since the
    # column widths are settled here
    texts, size = fit_titles([host for host, _, _ in hosts],
                             [width * w / sum(widths) - COLUMN_PADDING for w in widths])
    for title, text in zip(figure.layout.annotations, texts):
        title.update(text=text, font=dict(size=size))

    return figure, hosts


def add_band(figure, hosts, field, palette, labelled, row=2, legend='legend2',
             unknown=UNKNOWN_COLOR):
    '''
    A strip under each column saying one thing about each parasite of it -- the
    taxonomic group it belongs to, the niche it occupies -- one segment per parasite in
    the colour the palette gives that value.
    '''
    for column, (host, host_df, names) in enumerate(hosts, start=1):
        segments = dict(tuple(host_df.drop_duplicates('name').groupby(field, observed=True)))
        # the palette declares the order the values are read in; anything it does not name
        # follows behind
        ordered = ([v for v in palette if v in segments]
                   + [v for v in segments if v not in palette])
        for value in ordered:
            rows = segments[value]
            figure.add_trace(
                go.Bar(x=rows['name'], y=[1] * len(rows), name=value, width=1,
                       marker_color=palette.get(value, unknown),
                       # one slot per parasite: the values are a trace each for the legend's
                       # sake, and would otherwise be dealt half a column each
                       offsetgroup=field,
                       legend=legend, legendgroup=value,
                       showlegend=value not in labelled,
                       hovertemplate='%{x}' f'<extra>{value}</extra>'),
                row=row, col=column)
            labelled.add(value)
        figure.update_xaxes(categoryorder='array', categoryarray=names, row=row, col=column)
    figure.update_yaxes(visible=False, range=[0, 1], row=row)


def stack_bands(figure, upper=2, lower=3):
    '''
    Sit the two strips on one another. make_subplots leaves the same gap between every
    pair of rows, and the gap that keeps the strips apart from the columns is a line
    between the strips themselves, which say one thing each about the same parasite and
    read as one block under it.
    '''
    top = figure.get_subplot(lower, 1).yaxis.domain[1]
    bottom, ceiling = figure.get_subplot(upper, 1).yaxis.domain
    figure.update_yaxes(domain=[top, top + ceiling - bottom], row=upper)
    figure.update_yaxes(domain=[top + ceiling - bottom + BAND_GAP,
                                figure.get_subplot(1, 1).yaxis.domain[1]], row=1)


def add_group_and_niche_bands(figure, hosts, palette):
    '''
    Two strips under each column for the figure whose bars are split by localization:
    the taxonomic group of each parasite, and under it its niche.
    '''
    add_band(figure, hosts, 'group', palette, set(), row=2, legend='legend2')
    add_band(figure, hosts, 'niche', web_utils.NICHE_COLORS, set(), row=3, legend='legend3',
             unknown=web_utils.NICHE_COLORS[web_utils.UNKNOWN_NICHE])
    stack_bands(figure)
    figure.update_layout(height=520, margin=dict(t=175),
                         legend=dict(y=1.42, title_text='DeepLoc',
                                     title_font=dict(size=11)),
                         legend2=dict(orientation='h', yanchor='bottom', y=1.28, x=0,
                                      title_text='taxonomic group',
                                      title_font=dict(size=11), font=dict(size=11)),
                         legend3=dict(orientation='h', yanchor='bottom', y=1.14, x=0,
                                      title_text=web_utils.NICHE_TITLE,
                                      title_font=dict(size=11), font=dict(size=11)))
    # the names belong under the strips
    figure.update_xaxes(showticklabels=False, row=1)
    figure.update_xaxes(showticklabels=False, row=2)
    figure.update_xaxes(automargin=True, row=3)
    # the strips have pushed the counts off the figure
    figure.update_yaxes(automargin=True, row=1, col=1)


def add_niche_band(figure, hosts, labelled):
    '''
    The strip of niche under each column -- whether the parasite sits inside a host cell
    or outside it -- and the two legends a figure needs once it carries one, the colours
    of the bars above and the colours of the strip below being two keys to two different
    parts of it.
    '''
    add_band(figure, hosts, 'niche', web_utils.NICHE_COLORS, labelled,
             unknown=web_utils.NICHE_COLORS[web_utils.UNKNOWN_NICHE])
    figure.update_layout(margin=dict(t=125),
                         legend=dict(y=1.28, title_text='taxonomic group',
                                     title_font=dict(size=11)),
                         legend2=dict(orientation='h', yanchor='bottom', y=1.14, x=0,
                                      title_text=web_utils.NICHE_TITLE,
                                      title_font=dict(size=11), font=dict(size=11)))
    # the names belong under the strip, and the room they need is taken off the figure
    # rather than the margin
    figure.update_xaxes(showticklabels=False, row=1)
    figure.update_xaxes(automargin=True, row=2)


def style_host_columns(figure, y_title):
    '''
    The layout the two figures share: the parasite names under each column, the quantity
    named once down the left, and the taxonomic groups as the legend.
    '''
    figure.update_layout(height=470, plot_bgcolor='white',
                         margin=dict(l=0, r=0, t=95, b=10),
                         legend=dict(orientation='h', yanchor='bottom', y=1.14, x=0,
                                     title_text='', font=dict(size=11)))
    figure.update_xaxes(tickangle=-60, showgrid=False, tickfont=dict(size=11))
    figure.update_yaxes(showgrid=True, gridcolor='#f0f0f0', zerolinecolor='#e0e0e0')
    figure.update_yaxes(title_text=y_title, row=1, col=1)

    return figure


def classify_side(df, side):
    '''
    Which localization class DeepLoc puts the protein of `side` in, read on the classes
    that side was filtered on.
    '''
    if side == 'source':
        return web_utils.classify_localisation(df, web_utils.SURFACE_CLASSES,
                                               web_utils.BOTH_SURFACE)

    surface = pd.Series(web_utils.NOT_SURFACE, index=df.index)
    for niche, rows in df.groupby('niche'):
        surface.loc[rows.index] = web_utils.classify_localisation(
            rows, web_utils.niche_classes(niche), web_utils.SEVERAL)

    return surface


@st.cache_data(show_spinner=False)
def get_interaction_localisations(df, data_dir, side):
    '''
    The overview predictions with the localization class of one of the two proteins each
    interaction is between, so the bars counting the interactions can be split by where
    that protein sits.
    '''
    localisations = web_utils.load_deeploc_localisations(data_dir)
    if localisations.empty:
        return None

    df = df.merge(localisations, left_on=side, right_on='protein', how='left')
    df['surface'] = classify_side(df, side)

    return df


@st.cache_data(show_spinner=False)
def get_interactor_proteins(data_dir, config, side):
    '''
    One side of the predicted interactions, protein by protein, with what DeepLoc says
    about where each protein sits: the probability of each localization class the side
    was filtered on, and which of them the protein was called for -- one class, several
    of them, or none.
    '''
    localisations = web_utils.load_deeploc_localisations(data_dir)
    if localisations.empty:
        return None

    predictions = web_utils.load_predictions(data_dir)
    frames = []
    for taxid, host in config['hosts'].items():
        frame = predictions.loc[predictions['taxid2'] == str(taxid),
                                ['taxid1_label', side]].drop_duplicates()
        if not frame.empty:
            frames.append(frame.assign(host=host['label']))
    df = pd.concat(frames, ignore_index=True).merge(localisations, left_on=side,
                                                    right_on='protein', how='inner')

    groups = {p['label']: p.get('group', UNKNOWN_GROUP) for p in config['parasites'].values()}
    order = {g: i for i, g in enumerate(config.get('parasite_groups', {}))}
    df['group'] = df['taxid1_label'].map(lambda p: groups.get(p, UNKNOWN_GROUP))
    df['group_rank'] = df['group'].map(lambda g: order.get(g, len(order)))
    df['name'] = df['taxid1_label'].map(short_name)
    df['niche'] = df['taxid1_label'].map(web_utils.get_niches(config)).fillna(
        web_utils.UNKNOWN_NICHE)

    df['surface'] = classify_side(df, side)

    return df


@st.cache_data(show_spinner=False)
def generate_surface_scores_per_parasite(proteins, palette, width, score, cutoff, y_title,
                                         point_size=2.5):
    '''
    The spread of the probability itself, before it is a class: one box per parasite
    over the proteins of one surface class, in the same columns and colours as the
    figures above it.
    '''
    figure, hosts = host_columns(proteins, width, bands=1)
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
    figure.update_yaxes(range=[score_floor(cutoff, proteins[score].min()), 1], automargin=True,
                        row=1, col=1)
    # in the row of the columns rather than across the strip under them
    figure.add_hline(y=cutoff, line_width=1, line_dash='dot', line_color='#969696',
                     row=1, col='all')
    add_niche_band(figure, hosts, set())

    return figure


@st.cache_data(show_spinner=False)
def generate_host_score_boxes(proteins, width, point_size=3):
    '''
    How sure DeepLoc was of the host proteins, one column per localization class and one
    box per host inside it.
    '''
    # one row per protein and class it was called for; read from the probabilities, so a
    # protein over two cut-offs is a row under each
    columns = []
    for surface_class in web_utils.HOST_CLASSES:
        score_column = web_utils.DEEPLOC_SCORES[surface_class]
        if score_column not in proteins:
            continue
        # only the proteins some parasite was allowed to meet in that class
        allowed = proteins[proteins['niche'].map(
            lambda niche: surface_class in web_utils.niche_classes(niche))]
        counted = allowed.drop_duplicates(['host', 'protein'])
        in_class = counted[counted[score_column] > web_utils.DEEPLOC_CUTOFFS[surface_class]]
        columns.append(in_class.assign(surface=surface_class, score=in_class[score_column]))
    scored = pd.concat(columns, ignore_index=True)
    # here a column is a localization class and what is in it are the hosts, so the two are
    # swapped
    host_order = {host: rank for rank, host in enumerate(scored['host'].unique())}
    # four classes to a figure leave no room for full names
    scored['name'] = scored['host'].map(short_name)
    scored['taxid1_label'] = scored['host']
    scored['group_rank'] = scored['host'].map(host_order)
    scored['host'] = scored['surface']

    figure, classes = host_columns(scored, width)
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
        # in the column it applies to, not across the figure
        figure.add_hline(y=web_utils.DEEPLOC_CUTOFFS[surface_class], line_width=1,
                         line_dash='dot', line_color=SURFACE_LINE_COLORS[surface_class],
                         row=1, col=column)

    figure = style_host_columns(figure, 'P(the assigned localization)')
    # the scale clears the cut-offs of the classes actually drawn
    drawn = [web_utils.DEEPLOC_CUTOFFS[c] for c in scored['surface'].unique()]
    figure.update_yaxes(range=[score_floor(*drawn, scored['score'].min()), 1],
                        automargin=True, row=1, col=1)
    figure.update_xaxes(tickangle=0, automargin=True)

    return figure


@st.cache_data(show_spinner=False)
def generate_interactions_per_parasite(df, palette, width, score, classes=None):
    '''
    How many interactions are predicted for each parasite at or above a confidence, in
    one column per host, split by where DeepLoc puts the protein on one side of each of
    them.
    '''
    if classes is None:
        return generate_interactions_by_group(df, palette, width, score)

    # a class nothing was called for is still a legend entry; NOT_SURFACE is drawn only
    # where there is one
    classes = list(classes)
    if (df['surface'] == web_utils.NOT_SURFACE).any():
        classes.append(web_utils.NOT_SURFACE)

    figure, hosts = host_columns(df, width, bands=2)
    labelled = set()
    for column, (host, host_df, names) in enumerate(hosts, start=1):
        kept = host_df[host_df['weight'] >= score]
        counts = (kept.pivot_table(index='name', columns='surface', values='weight',
                                   aggfunc='count', fill_value=0)
                  .reindex(index=names, columns=classes, fill_value=0).fillna(0))
        # a parasite the threshold emptied has no column and is drawn as the absence it is
        total = counts.sum(axis=1).replace(0, pd.NA)
        for surface_class in classes:
            figure.add_trace(
                go.Bar(x=counts.index, y=counts[surface_class],
                       name=surface_class, marker_color=SURFACE_COLORS[surface_class],
                       customdata=pd.DataFrame({'n': counts[surface_class],
                                                'share': counts[surface_class] / total}),
                       legendgroup=surface_class, showlegend=surface_class not in labelled,
                       hovertemplate='%{x}<br>%{customdata[0]} predicted interactions '
                                     '(%{customdata[1]:.0%} of the parasite)'
                                     f'<extra>{surface_class}</extra>'),
                row=1, col=column)
            labelled.add(surface_class)
        figure.update_xaxes(categoryorder='array', categoryarray=names, row=1, col=column)

    figure = style_host_columns(figure, 'predicted interactions')
    figure.update_layout(barmode='stack', bargap=0.2)
    add_group_and_niche_bands(figure, hosts, palette)

    return figure


def generate_interactions_by_group(df, palette, width, score):
    '''
    The same counts with nothing to split them by: whole bars in the colour of the
    parasite's taxonomic group, with the niche in the one strip under them.
    '''
    figure, hosts = host_columns(df, width, bands=1)
    labelled_groups = set()
    for column, (host, host_df, names) in enumerate(hosts, start=1):
        kept = host_df[host_df['weight'] >= score]
        counts = kept.groupby(['name', 'group'], observed=True).size().reset_index(name='count')
        for group, rows in counts.groupby('group', observed=True):
            figure.add_trace(
                go.Bar(x=rows['name'], y=rows['count'], name=group,
                       marker_color=palette.get(group, UNKNOWN_COLOR),
                       # a bar per parasite, so the bars stand over the segments of the
                       # strip below them
                       offsetgroup='parasite',
                       legendgroup=group, showlegend=group not in labelled_groups,
                       hovertemplate='%{x}<br>%{y} predicted interactions'
                                     f'<extra>{group}</extra>'),
                row=1, col=column)
            labelled_groups.add(group)
        figure.update_xaxes(categoryorder='array', categoryarray=names, row=1, col=column)

    figure = style_host_columns(figure, 'predicted interactions')
    add_niche_band(figure, hosts, set())

    return figure


@st.cache_data(show_spinner=False)
def generate_confidence_per_parasite(df, palette, width, score):
    '''
    The spread of the confidence score of each parasite's predicted interactions, in the
    same columns and colours as the counts above, so the two figures are read together:
    a parasite with many interactions and a low box has many weakly supported ones.
    '''
    figure, hosts = host_columns(df, width, bands=1)
    labelled = set()
    for column, (host, host_df, names) in enumerate(hosts, start=1):
        for group, rows in host_df.groupby('group', observed=True):
            # the box statistics rather than every score, which would send each prediction
            # to the browser; hazen is the quartile plotly computes, and without points the
            # whiskers run to the extremes
            stats = rows.groupby('name', observed=True)['weight'].apply(
                lambda w: pd.Series(np.quantile(w, BOX_QUANTILES, method='hazen'),
                                    index=BOX_QUANTILES)).unstack()
            figure.add_trace(
                go.Box(x=stats.index, lowerfence=stats[0], q1=stats[0.25],
                       median=stats[0.5], q3=stats[0.75], upperfence=stats[1], name=group,
                       marker_color=palette.get(group, UNKNOWN_COLOR),
                       line_width=1, boxpoints=False,
                       legendgroup=group, showlegend=group not in labelled,
                       hovertemplate='%{x}<br>score %{y}' f'<extra>{host}</extra>'),
                row=1, col=column)
            labelled.add(group)
        figure.update_xaxes(categoryorder='array', categoryarray=names, row=1, col=column)

    figure = style_host_columns(figure, 'confidence score')
    # in the row of the columns rather than across the strip under them
    figure.add_hline(y=score, line_width=1, line_dash='dot', line_color='#969696',
                     row=1, col='all')
    # the line is drawn in every column and named in the first
    figure.add_annotation(text='counted above', x=0, y=score, xref='x domain', yref='y',
                          xanchor='left', yanchor='bottom', showarrow=False,
                          font=dict(size=10, color='#969696'), row=1, col=1)
    add_niche_band(figure, hosts, set())

    return figure


def og_label(names):
    '''
    The name a host orthologous group is drawn under: the gene symbols of its proteins,
    since the group id names nothing to read. The symbols are upper-cased so the same
    group reads the same way whichever hosts it is reached in.
    '''
    return web_utils.wrap_label(', '.join(sorted(set(str(n).upper() for n in names))),
                                LABEL_CHARS)


@st.cache_data(show_spinner=False)
def get_top_shared_ogs(df, score, top=TOP_OGS):
    '''
    The host orthologous groups the most parasites are predicted to interact with, across
    every host: the orthologous group is what the same protein of two hosts has in
    common. One row per parasite, host and orthologous group.
    '''
    kept = df[df['weight'] >= score]
    dots = kept.groupby(['host', 'taxid1_label', 'group2'], observed=True).agg(
        degree=('source', 'nunique'),
        # a tuple: streamlit hashes a dataframe through pandas, which cannot factorize a
        # list
        targets=('target_name', lambda n: tuple(sorted(set(str(x) for x in n))))).reset_index()
    parasites = dots.groupby('group2')['taxid1_label'].nunique()
    pairs = dots.groupby('group2').size()
    ranked = pd.DataFrame({'parasites': parasites, 'pairs': pairs})
    ranked = ranked[ranked['parasites'] > 1].sort_values(['parasites', 'pairs'],
                                                         ascending=False, kind='stable')
    if ranked.empty:
        return None

    ogs = list(ranked.head(top).index)
    dots = dots[dots['group2'].isin(ogs)].copy()
    dots['parasites'] = dots['group2'].map(ranked['parasites'])
    dots['pairs'] = dots['group2'].map(ranked['pairs'])
    dots['host proteins'] = dots['targets'].map(', '.join)
    # named by every symbol the group carries in any host, not only the ones this dot
    # reaches
    labels = {og: og_label(n for names in rows['targets'] for n in names)
              for og, rows in dots.groupby('group2')}
    dots['og'] = dots['group2'].map(labels)
    # the columns of the skeleton: the same parasite annotations the other figures carry
    parasite_rows = kept[['host', 'taxid1_label', 'name', 'group', 'group_rank',
                          'niche']].drop_duplicates()
    dots = dots.merge(parasite_rows, on=['host', 'taxid1_label'])

    return dots, [labels[f] for f in ogs], len(ranked)


@st.cache_data(show_spinner=False)
def get_shared_og_table(df, score):
    '''
    Every host orthologous group more than one parasite is predicted to interact with,
    for download, since the dot plot shows only the top of it. One row per group.
    '''
    kept = df[df['weight'] >= score]
    pairs = kept.groupby(['group2', 'host', 'taxid1_label'], observed=True)['target_name'].agg(
        lambda n: ', '.join(sorted(set(str(x) for x in n)))).reset_index()
    rows = []
    for og, reached in pairs.groupby('group2', observed=True):
        parasites = reached['taxid1_label'].nunique()
        if parasites < 2:
            continue
        symbols = sorted(set(s.upper() for names in reached['target_name']
                             for s in names.split(', ')))
        rows.append({
            'Orthologous group': og,
            'Gene symbols': ', '.join(symbols),
            'Parasites': parasites,
            'Host-parasite pairs': len(reached),
            'Interactions': '; '.join(f'{p} in {h}: {n}' for h, p, n
                                      in reached[['host', 'taxid1_label',
                                                  'target_name']].values)})

    return pd.DataFrame(rows).sort_values(['Parasites', 'Host-parasite pairs'],
                                          ascending=False, kind='stable')


@st.cache_data(show_spinner=False)
def generate_shared_og_dots(df, dots, ogs, palette, width):
    '''
    A dot wherever a parasite is predicted to interact with a protein of one of the
    orthologous groups, in the columns the figures above use, the groups ordered by how
    many parasites reach them so the most widely reached is the top row.
    '''
    # a wrapped name opens its row up
    ys = web_utils.row_positions(ogs, DOT_ROW)
    height = round(DOT_ROW * (ys[-1] + 1)) + DOT_CHROME
    # the strips keep the height they have under the other figures
    figure, hosts = host_columns(df, width, bands=2,
                                 band_height=BAND_HEIGHT * 470 / height)
    rows_of = dict(zip(ogs, ys))
    dots = dots.assign(y=dots['og'].map(rows_of))
    # the largest dot is as wide as a column has room for, its area standing for the
    # largest degree
    columns = sum(max(len(names), MIN_COLUMN) for _, _, names in hosts)
    room = width - (LABEL_CHAR * web_utils.label_width(ogs) + LABEL_PADDING)
    size = min(15, max(5, room / columns))
    sizeref = max(1, dots['degree'].max()) / size ** 2

    labelled = set()
    for column, (host, host_df, names) in enumerate(hosts, start=1):
        host_dots = dots[dots['host'] == host]
        # one trace per group, so the groups are the legend
        for group in [g for g in list(palette) + [UNKNOWN_GROUP]
                      if g in set(host_dots['group'])]:
            rows = host_dots[host_dots['group'] == group]
            figure.add_trace(
                go.Scatter(x=rows['name'], y=rows['y'], mode='markers', name=group,
                           marker=dict(color=palette.get(group, UNKNOWN_COLOR),
                                       size=rows['degree'], sizemode='area',
                                       sizeref=sizeref, sizemin=4, line=dict(width=0)),
                           legendgroup=group, showlegend=group not in labelled,
                           customdata=rows[['og', 'parasites', 'pairs', 'degree',
                                            'host proteins', 'group2']].to_numpy(),
                           hovertemplate='%{customdata[0]} (%{customdata[5]})<br>'
                                         'parasites reaching it: %{customdata[1]}, in '
                                         '%{customdata[2]} host-parasite pairs<br>'
                                         'proteins of %{x} reaching it: %{customdata[3]}<br>'
                                         'host proteins reached: %{customdata[4]}'
                                         '<extra></extra>'),
                row=1, col=column)
            labelled.add(group)
        # the range is set rather than left to the markers, which plotly pads at either
        # end, so the columns stand over the segments of the strip below them
        figure.update_xaxes(categoryorder='array', categoryarray=names,
                            range=[-0.5, len(names) - 0.5], row=1, col=column)

    figure = style_host_columns(figure, 'host orthologous group')
    # the dots are the legend of the groups already, so the group strip adds no entry;
    # the two share a legendgroup, and toggle together
    add_band(figure, hosts, 'group', palette, set(palette) | {UNKNOWN_GROUP},
             row=2, legend='legend2')
    add_band(figure, hosts, 'niche', web_utils.NICHE_COLORS, set(), row=3,
             legend='legend3', unknown=web_utils.NICHE_COLORS[web_utils.UNKNOWN_NICHE])
    stack_bands(figure)
    # reversed, so the most-shared group is the top row; the ticks are set on every
    # column so the grids line up
    figure.update_yaxes(range=[ys[-1] + 0.5, -0.5], tickmode='array',
                        tickvals=ys, ticktext=ogs,
                        zeroline=False, row=1)
    # the names are longer than the counts the margin was set for
    figure.update_yaxes(automargin=True, row=1, col=1)
    figure.update_xaxes(showgrid=True, gridcolor='#f0f0f0', row=1)
    # the names belong under the strips
    figure.update_xaxes(showticklabels=False, row=1)
    figure.update_xaxes(showticklabels=False, row=2)
    figure.update_xaxes(automargin=True, row=3)
    # the legends are placed from the top of the figure rather than the plot, which the
    # rows have made tall
    figure.update_layout(height=height, margin=dict(t=2 * LEGEND_ROW),
                         legend=dict(yref='container', yanchor='top', y=1,
                                     title_text='taxonomic group',
                                     title_font=dict(size=11)),
                         legend3=dict(orientation='h', x=0, yref='container',
                                      yanchor='top', y=1 - LEGEND_ROW / height,
                                      title_text=web_utils.NICHE_TITLE,
                                      title_font=dict(size=11), font=dict(size=11)))

    return figure


def common_name(host):
    '''`Homo sapiens (human)` as `human`.'''
    return host.split('(')[-1].rstrip(')') if '(' in host else host


@st.cache_data(show_spinner=False)
def get_host_group_presence(df, orthologs):
    '''
    The host orthologous groups each host has a protein in: those of the EggNOG members
    table, and every group a prediction reaches in that host, which the table may lack.
    Without the table every host is taken to have every group.
    '''
    if orthologs is None:
        everything = set(df['group2'].astype(str))
        return {host: everything for host in df['host'].unique()}

    present = {host: set(rows['group2'].astype(str))
               for host, rows in df.groupby('host', observed=True)}
    labels = {str(taxid): host['label'] for taxid, host in config['hosts'].items()}
    for taxid, rows in orthologs.groupby('taxid'):
        if str(taxid) in labels:
            present.setdefault(labels[str(taxid)], set()).update(rows['group'])

    return present


def pair_label(host, parasite):
    '''A parasite in one of its hosts, as a row of the matrix is named.'''
    return f'{parasite} ({common_name(host)})'


@st.cache_data(show_spinner=False)
def get_cross_host_similarity(df, score, hosts, present, order_by):
    '''
    How alike the host orthologous groups reached by each pair of host-parasite pairs are,
    as the Jaccard similarity of the two sets. Only the groups both hosts have a protein in
    are counted, so a group one host lacks is not read as a difference between the
    parasites. One row per parasite and host, in the order OVERLAP_ORDERS gives
    `order_by`.
    '''
    kept = df[(df['weight'] >= score) & df['host'].isin(hosts)]
    reached = {key: set(rows['group2'].astype(str))
               for key, rows in kept.groupby(['host', 'taxid1_label'], observed=True)}
    if len(reached) < 3:
        return None

    host_rank = {h: i for i, h in enumerate(hosts)}
    rows = kept.drop_duplicates(['host', 'taxid1_label'])
    rows = rows.assign(niche_rank=rows['niche'].map(web_utils.niche_rank),
                       host_rank=rows['host'].map(host_rank)).sort_values(
        OVERLAP_ORDERS[order_by], kind='stable')
    keys = list(zip(rows['host'], rows['taxid1_label']))

    shared = np.full((len(keys), len(keys)), np.nan)
    similarity = np.full((len(keys), len(keys)), np.nan)
    for i, (host_a, parasite_a) in enumerate(keys):
        for j, (host_b, parasite_b) in enumerate(keys[:i]):
            common = present.get(host_a, set()) & present.get(host_b, set())
            first = reached[host_a, parasite_a] & common
            second = reached[host_b, parasite_b] & common
            union = len(first | second)
            shared[i, j] = shared[j, i] = len(first & second)
            similarity[i, j] = similarity[j, i] = shared[i, j] / union if union else 0

    labels = [pair_label(h, p) for h, p in keys]

    return (pd.DataFrame(similarity, index=labels, columns=labels),
            pd.DataFrame(shared, index=labels, columns=labels),
            list(rows['group']), list(rows['niche']),
            [common_name(h) for h in rows['host']], keys)


def og_table(edges, groups, keys):
    '''
    One row per orthologous group, named by the gene symbols the two parasites reach in
    it, with a column per host-parasite pair counting the parasite proteins reaching it.
    '''
    edges = edges[edges['group2'].astype(str).isin(groups)]
    table = pd.DataFrame({'Orthologous group': sorted(groups)})
    symbols = edges.groupby(edges['group2'].astype(str))['target_name'].agg(
        lambda names: ', '.join(sorted(set(str(n).upper() for n in names))))
    table['Gene symbols'] = table['Orthologous group'].map(symbols).fillna('')
    columns = []
    for host, parasite in keys:
        rows = edges[(edges['host'] == host) & (edges['taxid1_label'] == parasite)]
        degree = rows.groupby(rows['group2'].astype(str))['source'].nunique()
        column = f'{short_name(parasite)} ({common_name(host)}) proteins'
        table[column] = table['Orthologous group'].map(degree).fillna(0).astype(int)
        columns.append(column)

    # the groups reached by the most parasite proteins first
    return table.assign(_reach=table[columns].sum(axis=1)).sort_values(
        ['_reach', 'Orthologous group'], ascending=[False, True],
        kind='stable').drop(columns='_reach')


@st.dialog('Host orthologous groups of a pair of parasites', width='large')
def show_cross_host_dialog(first, second, df, score, present):
    '''
    The orthologous groups behind one cell of the cross-host matrix, split by why a group
    is or is not shared: reached by both, reached by one though both hosts have it, or
    reached by one in a host the other has no protein of the group in.
    '''
    keys = [first, second]
    edges = df[df['weight'] >= score]
    edges = edges[pd.Series(list(zip(edges['host'], edges['taxid1_label'])),
                            index=edges.index).isin(keys)]
    reached = {key: set(rows['group2'].astype(str))
               for key, rows in edges.groupby(['host', 'taxid1_label'], observed=True)}
    first_groups, second_groups = reached.get(first, set()), reached.get(second, set())
    common = present.get(first[0], set()) & present.get(second[0], set())

    shared = first_groups & second_groups & common
    one_parasite = (first_groups ^ second_groups) & common
    no_ortholog = (first_groups | second_groups) - common

    st.markdown(f'**{pair_label(*first)}** and **{pair_label(*second)}**')
    st.caption('Host orthologous groups the two parasites are predicted to interact with, '
               'at the confidence set on the page. A column counts the proteins of that '
               'parasite predicted to reach the group. The matrix counts only the first two '
               'tabs: a group one host has no protein in cannot be reached there.')
    tabs = st.tabs([f'Shared ({len(shared)})',
                    f'Reached by one parasite ({len(one_parasite)})',
                    f'No ortholog in the other host ({len(no_ortholog)})'])
    for tab, groups, empty in [
            (tabs[0], shared, 'The two parasites share no host orthologous group.'),
            (tabs[1], one_parasite, 'Every group both hosts have is reached by both '
                                    'parasites or by neither.'),
            (tabs[2], no_ortholog, 'Every group either parasite reaches has a protein in '
                                   'both hosts.')]:
        with tab:
            if groups:
                st.dataframe(og_table(edges, groups, keys), width='stretch', height=400,
                             hide_index=True)
            else:
                st.info(empty)


st.caption('Host-parasite protein-protein interactions predicted by '
           'orthology transfer. This page presents every prediction per host, '
           '**Host-parasite network** shows the network of one host-parasite pair, '
           '**Parasites of a host** compares the parasites of a single host, and **Multi-host '
           'parasites** compares the hosts of multi-host parasites.')
st.markdown("---")

st.subheader('Hosts and parasites')
overview = get_overview_predictions(data_dir, config)
hosts_parasites.show(config)
st.markdown("---")

# measured once here so the column titles can be fitted; page_width answers with a default
# until the browser has replied
page = web_utils.page_width()
parasite_palette = config.get('parasite_groups', {})

# one slider for the page, in the middle column so it does not take the whole width
with st.columns(3)[1]:
    score = web_utils.score_slider(
        help='Interactions predicted below this confidence are left out '
             'of the counts and of the localization split drawn on them. The '
             'boxplots keep every prediction: the confidence figure draws '
             'this threshold as a line instead, and the localization ones '
             'stand on too few proteins to be thresholded as well.')
    # what the threshold keeps, under the slider
    kept = int((overview['weight'] >= score).sum())
    st.caption(f'{kept:,} of {len(overview):,} interactions ({kept / len(overview):.0%}) '
               f'at or above {score:g}.')

# the description in the heading depends on the toggles but sits above them
heading = st.empty()
st.caption('How many interactions are predicted for each parasite, in each of its hosts.')
# narrow columns: a segmented control left to itself is stretched over the page
with st.columns([1.4, 1, 1.6])[0]:
    split_side = st.segmented_control(
        'Colour by', list(SPLIT_SIDES), default=list(SPLIT_SIDES)[0],
        key='split_side',
        help='The taxonomic group of the parasite, or the localization DeepLoc 2 assigns '
             'the host or the parasite protein of each interaction, which splits the bars.')
side, split_classes = SPLIT_SIDES[split_side or list(SPLIT_SIDES)[0]]
interactions = None if side is None else get_interaction_localisations(overview, data_dir,
                                                                       side)

if interactions is None:
    description = ('Predicted interactions per parasite at or above the confidence set above, '
                   'grouped by host and coloured by parasite taxonomic group. ' + NICHE_STRIP)
elif side == 'target':
    description = ('Predicted interactions per parasite at or above the confidence set above, '
                   'grouped by host and split by the localization DeepLoc 2 '
                   'assigns the host protein of each interaction. ' + BANDS_STRIP)
else:
    description = ('Predicted interactions per parasite at or above the confidence set above, '
                   'grouped by host and split by the localization DeepLoc 2 assigns the '
                   'parasite protein of each interaction. ' + BANDS_STRIP)
heading.subheader("Number of predicted interactions per parasite", help=description)

st.plotly_chart(
    generate_interactions_per_parasite(overview if interactions is None else interactions,
                                       parasite_palette, page, score,
                                       classes=None if interactions is None else split_classes),
    width='stretch')

og_help = ('Each group is the EggNOG orthologous group of the host protein, which can hold '
           'several paralogous genes, and is named by the gene symbols of its proteins. A dot '
           'wherever a parasite is predicted to interact with a protein of the group, '
           "sized by the number of that parasite's proteins reaching it; hover a dot for the "
           'host proteins behind it. ' + BANDS_STRIP)
shared_ogs = get_top_shared_ogs(overview, score)
if shared_ogs is None:
    st.subheader("Most common host orthologous groups", help=og_help)
    st.info('No host orthologous group is reached by more than one parasite at this confidence.')
else:
    og_dots, ogs, shareable = shared_ogs
    shown = len(ogs)
    if shareable > shown:
        og_help = (f'The {shown} groups reached by the most parasites are shown, of the '
                   f'{shareable} reached by more than one at or above the confidence set '
                   'above. ' + og_help)
        selection = f'The {shown} host orthologous groups targeted by the most parasites.'
    else:
        og_help = ('Every group reached by more than one parasite at or above the confidence '
                   'set above is shown. ' + og_help)
        selection = f'The {shown} host orthologous groups targeted by more than one parasite.'
    st.subheader("Most common host orthologous groups", help=og_help)
    st.caption(selection)
    st.plotly_chart(generate_shared_og_dots(overview, og_dots, ogs,
                                            parasite_palette, page),
                    width='stretch')
    st.download_button(f'Download all {shareable} shared orthologous groups',
                       data=utils.convert_df(get_shared_og_table(overview, score)),
                       file_name='shared_host_orthologous_groups.tsv', mime='text/csv')

st.subheader(
    "Overlap of the host orthologous groups across hosts",
    help='How alike the host interactors of each pair of parasites are, across hosts as '
         'well as within one: the Jaccard similarity of the host orthologous groups the two '
         'reach, the groups both reach over the groups either reaches. A protein of one '
         'host is never the protein of another, so the orthologous group stands in for it. '
         'Only the groups both hosts have a protein in are counted, so a group missing '
         'from one host is not read as a difference between the parasites; the cell '
         'dialog lists those groups apart. A parasite predicted in several hosts has a row '
         'in each. Ordered by niche or by taxonomic group, a parasite or genus lines up '
         'across its hosts; ordered by host, each host is a block and the cells between '
         'two blocks compare the parasites of the two. Three strips run along each axis: the '
         'taxonomic group of the parasite, whether it lives inside a host cell or outside '
         'one, and the host. Interactions below the confidence set above are left out. '
         'Hover a cell for the count behind its ratio, or click it to see the groups.')
st.caption('How alike the host orthologous groups targeted by each pair of parasites are, '
           'in any of their hosts. Click a cell to see the groups.')
host_labels = list(dict.fromkeys(overview['host']))
hosts_column, order_column, _ = st.columns(3)
with hosts_column:
    overlap_hosts = st.multiselect(
        'Hosts', host_labels, default=host_labels, format_func=common_name,
        key='overlap_hosts',
        help='The hosts whose parasites the matrix compares, for instance human and mouse '
             'alone to set the Plasmodium species of the two side by side.')
with order_column:
    overlap_order = st.radio('Order the parasites by', list(OVERLAP_ORDERS), horizontal=True,
                             key='overlap_order',
                             help='Which annotation comes out in blocks along the axes. '
                                  'The strips show all three either way.')
host_orthologs = web_utils.load_host_orthologs(data_dir)
if host_orthologs is None:
    st.caption('Run `python scripts/build_host_orthologs.py` to count only the groups both '
               'hosts have; until then every group reached is counted.')
group_presence = get_host_group_presence(overview, host_orthologs)
cross_host = get_cross_host_similarity(overview, score, tuple(overlap_hosts),
                                       group_presence, overlap_order)
if cross_host is None:
    st.info('Fewer than three parasites of the selected hosts reach a host orthologous '
            'group at this confidence.')
else:
    similarity, shared_counts, clades, row_niches, row_hosts, row_keys = cross_host
    figure, click_layer = similarity_heatmap.generate(
        similarity, shared_counts,
        [(clades, {**parasite_palette, UNKNOWN_GROUP: UNKNOWN_COLOR}),
         (row_niches, {n: web_utils.NICHE_COLORS[n]
                       for n in web_utils.NICHE_ORDER + [web_utils.UNKNOWN_NICHE]}),
         (row_hosts, {common_name(h): config['hosts'][t]['color']
                      for t in config['hosts'] for h in [config['hosts'][t]['label']]})],
        page,
        x_names=[f'{short_name(p)} ({common_name(h)})' for h, p in row_keys],
        count_label='Shared orthologous groups',
        # held to about a screen, so the whole matrix is seen at once; with every host the
        # names go down to eight points, upright so the narrow columns still part them
        max_side=OVERLAP_SIDE, smallest_label=8, x_angle=-90)
    # remounted after every dialog: Streamlit drops an identical selection and plotly
    # reads a second click as a deselection
    nonce = st.session_state.get(CELL_NONCE_KEY, 0)
    clicked = st.plotly_chart(
        figure, width='content', on_select='rerun', selection_mode='points',
        key=f'cross_host_cells_{"_".join(overlap_hosts)}_{overlap_order}_{score}_{nonce}')
    cell = next((point for point in (clicked or {}).get('selection', {}).get('points', [])
                 if point.get('curve_number') == click_layer), None)
    if cell is not None:
        # the two host-parasite pairs are the row and column the marker sits on
        st.session_state[CELL_NONCE_KEY] = nonce + 1
        show_cross_host_dialog(row_keys[int(cell['y'])], row_keys[int(cell['x'])],
                               overview, score, group_presence)

st.subheader(
    "Interaction confidence scores",
    help='Boxplots of the distribution of confidence scores per parasite. Each score is '
         'the average of the experimental and database evidence channels of the KOG-KOG '
         'link the prediction was transferred from. Every prediction is counted here, '
         'whatever the slider is set to; the dotted line marks the threshold it selects '
         'instead. ' + NICHE_STRIP)
st.caption('How confident the predictions are for each parasite.')
st.plotly_chart(generate_confidence_per_parasite(overview, parasite_palette, page, score),
                width='stretch')

# one protein per row rather than one interaction, on every prediction and not at the
# threshold
host_proteins = get_interactor_proteins(data_dir, config, 'target')
parasite_proteins = get_interactor_proteins(data_dir, config, 'source')
# the membrane figure is drawn over the unicellular parasites alone
unicellular = None
if parasite_proteins is not None:
    unicellular = parasite_proteins[parasite_proteins['group'].isin(UNICELLULAR_GROUPS)]

if host_proteins is not None:
    st.subheader(
        "Localization probabilities of host proteins",
        help='Boxplots of the DeepLoc 2 probabilities of the host proteins for their '
             'assigned localization, one column per class and one box per host. The dotted '
             'line in each column marks the cut-off that class is called at, DeepLoc 2\'s own '
             'default threshold for the Accurate model: '
             + ', '.join(f'{web_utils.DEEPLOC_CUTOFFS[c]:.3f} for {c.lower()}'
                         for c in web_utils.HOST_CLASSES) + '. '
             'That is what the proteins were filtered on, so every point is above its own '
             'line. A protein over several cut-offs appears in each of those columns, in '
             'each at the probability of that class. The cytoplasm and nucleus columns hold '
             'the host proteins of the intracellular parasites, the only ones the filter '
             'allows them to. Each host protein is counted once per class, '
             'irrespective of the number of parasites reaching it, and every prediction is '
             'counted whatever the slider is set to.')
    st.caption('How confidently DeepLoc 2 places the host proteins where they were '
               'localized.')
    st.plotly_chart(generate_host_score_boxes(host_proteins, page), width='stretch')

if parasite_proteins is not None:
    st.subheader(
        "Localization probabilities of extracellular parasite proteins",
        help='Boxplots of the DeepLoc 2 extracellular probability for proteins assigned '
             'extracellular or both classes for each parasite, over every prediction whatever '
             'the slider is set to. The dotted line marks the cut-off, DeepLoc 2\'s default '
             'threshold for extracellular under the Accurate model '
             f'({web_utils.DEEPLOC_CUTOFFS[web_utils.EXTRACELLULAR]:.3f}), which is '
             'what the secretome filter kept these proteins on. Individual proteins are '
             'shown as points behind each box; for parasites with a hundred or more proteins '
             'the points are read as density. ' + NICHE_STRIP)
    st.caption('How confidently DeepLoc 2 calls the parasite proteins secreted.')
    st.plotly_chart(
        generate_surface_scores_per_parasite(
            parasite_proteins[parasite_proteins['surface'].isin(
                [web_utils.EXTRACELLULAR, web_utils.BOTH_SURFACE])],
            parasite_palette, page, 'extracellular',
            web_utils.DEEPLOC_CUTOFFS[web_utils.EXTRACELLULAR],
            # forty-five columns and up to 180 proteins in one, so the smallest dot that
            # still carries colour
            'P(extracellular)', point_size=2),
        width='stretch')

    st.subheader(
        "Localization probabilities of membrane parasite proteins",
        help='The equivalent for parasite proteins assigned to cell membrane or both '
             'classes, scored on that probability, again over every prediction. Only '
             'unicellular parasites are represented, the secretome filter admitting a '
             'multicellular parasite nothing but its secreted proteins. Individual '
             'proteins are shown as points; boxes over very few proteins (two for '
             '*V. corneae*) should not be read as distributions. The '
             'dotted line marks the cut-off, DeepLoc 2\'s default for cell membrane under '
             f'the Accurate model ({web_utils.DEEPLOC_CUTOFFS[web_utils.CELL_MEMBRANE]:.3f}). '
             + NICHE_STRIP)
    st.caption('How confidently DeepLoc 2 places the parasite proteins in the cell '
               'membrane (unicellular parasites only).')
    st.plotly_chart(
        generate_surface_scores_per_parasite(
            unicellular[unicellular['surface'].isin(
                [web_utils.CELL_MEMBRANE, web_utils.BOTH_SURFACE])],
            parasite_palette, page, 'cell_membrane',
            web_utils.DEEPLOC_CUTOFFS[web_utils.CELL_MEMBRANE],
            # fewer parasites and proteins than the figure above, so the dots can be told
            # apart
            'P(cell membrane)', point_size=3),
        width='stretch')

st.markdown("---")


# Footer
with st.container():
    web_utils.footer()
