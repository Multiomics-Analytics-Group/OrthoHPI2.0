import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils
import web_utils
import hosts_parasites
import streamlit as st
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
# the parasite groups the membrane figure is drawn over. A multicellular parasite reaches
# its host with secreted proteins alone -- the secretome filter keeps it nothing else -- so
# its split is the filter and not the parasite; only the unicellular groups have both
# surface classes open to them
UNICELLULAR_GROUPS = ('Apicomplexa', 'Kinetoplastida', 'Other protozoa', 'Microsporidia')
UNKNOWN_COLOR = '#999999'
# A host with two parasites still has to fit two labels under its column, so the columns
# are sized as if every host had at least this many parasites. Human has 35 against the
# two of pig, and strictly proportional columns leave pig a strip its labels overrun.
MIN_COLUMN = 4
# what the strips under the columns are saying, written once for the figures that
# carry them. `niche` is the key in config.yml, but the term the literature uses for the
# split is the two words themselves, so the caption spells it out instead
INTRACELLULAR_NOTE = ('A parasite with an intracellular stage in the host counts as '
                      'intracellular, since the reach of that stage is the wider of the two.')
NICHE_STRIP = ('The strip below the columns indicates whether the parasite lives inside a '
               'host cell or outside one. ' + INTRACELLULAR_NOTE)
# and the same for the counts figure, which carries the taxonomic group under the columns as
# well, its bars having spent their colour on the localization classes
BANDS_STRIP = ('The upper strip below the columns indicates the taxonomic group of the '
               'parasite, the lower one whether it lives inside a host cell or outside '
               'one. ' + INTRACELLULAR_NOTE)
# the confidence the counts of the page are drawn at, and the range the slider spans, the
# same default and range as the network page: a parasite counted here then agrees with the
# network the reader opens next instead of being several times larger than it
MIN_SCORE, MAX_SCORE, DEFAULT_SCORE = 0.35, 0.9, 0.35
# the title over a column is the name of the host, and a host with two parasites has a
# column narrower than its own name. The name is written smaller, and on two lines, until
# it fits: plotly draws a subplot title at a size of its own choosing and centred, so two
# of them over two narrow columns are drawn across each other.
# About 0.55 of the font size a character, the names being proportional text, so this is an
# average rather than a measurement of the string
TITLE_CHAR = 0.55
# and the pixels of a column its title is not written across: the gap between two columns,
# so that two titles that both fill their columns still have a space between them
COLUMN_PADDING = 12
# plotly's own size for a subplot title, which is the largest one is written at here
TITLE_SIZE = 16
# and the smallest: below this the name of a host is no longer read, so a column too narrow
# for it keeps this size and the title is drawn a little wider than the column
SMALLEST_TITLE = 10
# share of the height of a figure taken by the strip of taxonomic group under its columns.
# Enough to read as a band of colour, not enough to be read as a quantity of its own
BAND_HEIGHT = 0.06
# and the share of it left blank between the columns and that strip
BAND_GAP = 0.04
# the colour of each localization class, kept in web_utils because the dot plot of the
# "Parasites of a host" page keys the strip beside its rows with the same palette
SURFACE_COLORS = web_utils.LOCALISATION_COLORS
# and what to outline a box of that class in, where the class is drawn as a box rather than
# as a bar: the pale shade is a fill and an outline drawn in it on a white background is an
# outline the reader has to look for. No entry for the mixed classes, drawn as bars only
SURFACE_LINE_COLORS = {'Extracellular': '#3690c0', 'Cell membrane': '#045a8d',
                       'Cytoplasm': '#e6550d', 'Nucleus': '#a63603'}
# the classes each side of an interaction is split over, and so the segments of a bar and
# the entries of the legend. The host side carries the four classes its filter reads; the
# parasite side the surface pair, which is all the secretome filter selected on
HOST_SPLIT_CLASSES = web_utils.HOST_CLASSES + (web_utils.SEVERAL,)
PARASITE_SPLIT_CLASSES = web_utils.SURFACE_CLASSES + (web_utils.BOTH_SURFACE,)
# An interaction has a protein at either end of it, so the counts figure can be split by
# either: what the parasite reaches its host with, or what of the host it reaches. Each is
# read on the classes its own filter allowed, which is why the two are not the same list.
# The label is what the toggle over the figure offers, the host side first as the side the
# parasites can be compared on -- the parasite side is partly the secretome filter
SPLIT_SIDES = {'Host proteins': ('target', HOST_SPLIT_CLASSES),
               'Parasite proteins': ('source', PARASITE_SPLIT_CLASSES)}
# and whether a column is a count or its own 100%. The counts run from 39 interactions to
# two and a half thousand, so the split of the smaller parasites is only legible in shares,
# while the counts are what the figure is about
BAR_SCALES = ('Counts', 'Share')
# how far under the lowest thing a probability scale has to show -- its cut-off, or a point
# below it -- the scale starts. Enough that the line and the points sitting on it are not
# drawn against the axis itself
SCALE_MARGIN = 0.05


def score_floor(*values):
    '''
    Where a probability scale starts: a little under the lowest of what has to be visible on
    it. The figures below are read within a class rather than across the whole 0 to 1, so a
    scale that starts at 0 is a scale whose boxes are flattened into the top of the figure.

    Called with the cut-off and the smallest probability drawn. Everything this pipeline
    builds is above its cut-off, both sides having been filtered on it, so the cut-off is
    what the scale clears; the smallest probability is there for the snapshot data
    directories, built before these cut-offs, where a point can fall under the line and
    clipping to the line would hide it.

    :param values: the probabilities the scale has to leave room for
    :return: the foot of the scale, never below 0
    '''
    return max(0.0, min(values) - SCALE_MARGIN)


def short_name(parasite):
    '''`Plasmodium vivax` as `P. vivax`, the same abbreviation the circos labels its
    arcs with, so the same parasite reads the same way on every page.'''
    return f'{parasite[0]}. {parasite.split(" ")[1]}'


def fit_titles(hosts, rooms):
    '''
    The names of the hosts, written to fit the columns they head.

    Three forms, each narrower than the last: the name as it is given, the same over two
    lines with the common name under the species, and the species abbreviated the way the
    parasites of the columns are. The one used is whichever can be written largest, and not
    the first that fits at all: a narrower form is worth the abbreviation when it buys a
    size, and the widest form that fits is often the one that fits at the smallest size
    there is.

    One form and one size for the whole figure rather than the widest each column could
    take on its own. The names are a row of labels across the top of one figure, and a row
    in which one host is abbreviated and another spelt out, at two sizes and on a different
    number of lines each, reads as four things rather than as four hosts.

    :param list hosts: the labels of the hosts, `Rattus norvegicus (rat)`, in column order
    :param list rooms: pixels each of those columns has for its name
    :return: (the texts, in plotly markup and in the same order, the size to write them at)
    '''
    def forms(host):
        species, _, common = host.partition(' (')
        common = f'({common}' if common else ''
        if not common:
            return [[species]] * 3

        return [[f'{species} {common}'], [species, common], [short_name(species), common]]

    def size_of(lines, room):
        return min(TITLE_SIZE, int(room / (TITLE_CHAR * max(len(line) for line in lines))))

    written = [[forms(host)[form] for host in hosts] for form in range(3)]
    # the size a form can be written at is the size its narrowest column allows, and the
    # form chosen is the one that comes out largest -- the earliest of them where two do,
    # a name spelt out being worth more than the same name abbreviated at the same size
    sizes = [min(size_of(lines, room) for lines, room in zip(form, rooms))
             for form in written]
    best = sizes.index(max(sizes))

    return ['<br>'.join(lines) for lines in written[best]], max(sizes[best], SMALLEST_TITLE)


def host_coverage_caption(data_dir, config):
    '''A compact summary of the post-filter host pools behind cross-host comparisons.'''
    eligible = web_utils.load_eligible_proteins(data_dir)
    if eligible is None:
        return None

    parts = []
    for taxid, host in config['hosts'].items():
        count = (eligible['taxid'].astype(str) == str(taxid)).sum()
        parts.append(f"{host['label']}: {count:,} proteins (TISSUES >= "
                     f"{host['tissue_cutoff']:g})")
    return 'Eligible host proteins after tissue and DeepLoc filtering — ' + '; '.join(parts) + '.'


@st.cache_data(show_spinner=False)
def get_overview_predictions(data_dir, config):
    '''
    Every predicted interaction, labelled with the host species it was predicted against
    and with the taxonomic group of its parasite. A parasite infecting more than one host
    appears once per host, since the point of this page is comparing those.

    :param str data_dir: directory holding predictions.parquet
    :param dict config: parsed configuration
    :return: one row per predicted interaction, with host, parasite, group, niche, weight
             and the two proteins it is between
    '''
    predictions = web_utils.load_predictions(data_dir)
    groups = {p['label']: p.get('group', UNKNOWN_GROUP) for p in config['parasites'].values()}
    order = {g: i for i, g in enumerate(config.get('parasite_groups', {}))}

    frames = []
    for taxid, host in config['hosts'].items():
        frame = predictions.loc[predictions['taxid2'] == str(taxid),
                                ['taxid1_label', 'weight', 'source', 'target']]
        if not frame.empty:
            frames.append(frame.assign(host=host['label']))
    df = pd.concat(frames, ignore_index=True)

    df['group'] = df['taxid1_label'].map(lambda p: groups.get(p, UNKNOWN_GROUP))
    # what puts the parasites of a host in the order of the circos: taxonomic group as
    # the configuration declares it, then name, with an unclassified parasite last
    df['group_rank'] = df['group'].map(lambda g: order.get(g, len(order)))
    df['name'] = df['taxid1_label'].map(short_name)
    df['niche'] = df['taxid1_label'].map(web_utils.get_niches(config)).fillna(
        web_utils.UNKNOWN_NICHE)

    return df


def host_columns(df, width, bands=0):
    '''
    The skeleton the figures are drawn on: one column per host, as wide as the number of
    parasites infecting it, sharing a y axis so a bar or a box can be compared straight
    across the hosts rather than only within one.

    :param df: overview predictions, as get_overview_predictions builds them
    :param float width: pixels the figure is drawn across, which is what the names over the
                        columns are fitted to. A column is a share of it, and the share a
                        host with two parasites gets is narrower than its own name
    :param int bands: shallow rows to leave under the columns for the strips add_band
                      draws, one per fact about the parasites the bars themselves cannot
                      carry -- their taxonomic group, their niche -- drawn from row 2 down
    :return: (figure, [(host, its predictions, its parasites in axis order), ...])
    '''
    hosts = []
    for host in df['host'].unique():
        host_df = df[df['host'] == host]
        parasites = host_df[['taxid1_label', 'name', 'group_rank']].drop_duplicates()
        parasites = parasites.sort_values(by=['group_rank', 'taxid1_label'], kind='stable')
        hosts.append((host, host_df, parasites['name'].tolist()))

    widths = [max(len(names), MIN_COLUMN) for _, _, names in hosts]
    figure = make_subplots(rows=1 + bands, cols=len(hosts), shared_yaxes=True,
                           column_widths=[w / sum(widths) for w in widths],
                           # the titles fill the top row, which is the row of the figure
                           subplot_titles=[host for host, _, _ in hosts],
                           row_heights=([1 - bands * BAND_HEIGHT] + [BAND_HEIGHT] * bands
                                        if bands else None),
                           # enough of a gap that a strip is read as a second thing about
                           # the columns rather than as the foot of the columns themselves
                           vertical_spacing=BAND_GAP, horizontal_spacing=0.02)

    # the titles as make_subplots left them are the first annotations of the figure, one per
    # column and in the order the columns were given. Rewritten rather than passed in
    # already fitted, since the room a title has is the share of the width its column came
    # out with, which is settled here
    texts, size = fit_titles([host for host, _, _ in hosts],
                             [width * w / sum(widths) - COLUMN_PADDING for w in widths])
    for title, text in zip(figure.layout.annotations, texts):
        title.update(text=text, font=dict(size=size))

    return figure, hosts


def add_band(figure, hosts, field, palette, labelled, row=2, legend='legend2',
             unknown=UNKNOWN_COLOR):
    '''
    A strip under each column saying one thing about each parasite of it -- the taxonomic
    group it belongs to, the niche it occupies -- one segment per parasite in the colour
    the palette gives that value. The segments touch, so the strip reads as a strip and the
    parasites of a host sharing a value are a block rather than a row of bars.

    A strip gets a legend of its own, under the one naming the colours of the bars: the two
    say different things about different parts of the figure, and in a single row of entries
    the values of the strip read as more of what the bars are split into.

    :param figure: a figure host_columns built with a band to spare, modified in place
    :param list hosts: the (host, its rows, its parasites in axis order) of host_columns
    :param str field: the column of the frame the strip is drawn from
    :param dict palette: {value of that column: colour}
    :param set labelled: values already in the legend, added to as they are drawn
    :param int row: the row of the figure to draw the strip in, counting the columns as 1
    :param str legend: the legend its entries belong to
    :param str unknown: colour for a value the palette does not name
    '''
    for column, (host, host_df, names) in enumerate(hosts, start=1):
        segments = dict(tuple(host_df.groupby(field, observed=True)))
        # the palette declares the order the values are read in -- the clades in the order
        # of the circos, the niches from outside the host cell to inside it -- and anything
        # it does not name follows behind
        ordered = ([v for v in palette if v in segments]
                   + [v for v in segments if v not in palette])
        for value in ordered:
            rows = segments[value]
            figure.add_trace(
                go.Bar(x=rows['name'], y=[1] * len(rows), name=value, width=1,
                       marker_color=palette.get(value, unknown),
                       # one segment per parasite and one slot per parasite: the values are
                       # a trace each for the legend's sake, and traces left in slots of
                       # their own are dealt half a column each and drawn off centre
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
    Sit the two strips on one another. make_subplots leaves the same gap between every pair
    of rows, and the gap that keeps the strips apart from the columns is a line between the
    strips themselves, which say one thing each about the same parasite and read as one
    block under it. The upper strip drops onto the lower one and the columns take back the
    room it left, so the figure keeps its height and its gap under the columns.

    :param figure: a figure host_columns built with two bands, modified in place
    :param int upper: the row of the strip that moves down
    :param int lower: the row of the strip it comes to rest on
    '''
    top = figure.get_subplot(lower, 1).yaxis.domain[1]
    bottom, ceiling = figure.get_subplot(upper, 1).yaxis.domain
    figure.update_yaxes(domain=[top, top + ceiling - bottom], row=upper)
    figure.update_yaxes(domain=[top + ceiling - bottom + BAND_GAP,
                                figure.get_subplot(1, 1).yaxis.domain[1]], row=1)


def add_group_and_niche_bands(figure, hosts, palette):
    '''
    Two strips under each column for the figure whose bars are split by localization: the
    taxonomic group of each parasite, and under it its niche. The colour of the bars is
    spent on the classes, so both facts about the parasite itself move below the columns,
    the clades of a host reading as blocks of colour in the upper strip.

    Three legends then, one per key: the classes the bars are split into on top, the clades
    of the upper strip under it, and the niche of the lower strip last, in the order the
    three parts of the figure are met going down it.

    :param figure: a figure host_columns built with two bands to spare, modified in place
    :param list hosts: the (host, its rows, its parasites in axis order) of host_columns
    :param dict palette: {taxonomic group: colour}
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
    # the names belong under the strips, which are the foot of the figure now, and a column
    # labelled three times is a column labelled twice too often
    figure.update_xaxes(showticklabels=False, row=1)
    figure.update_xaxes(showticklabels=False, row=2)
    figure.update_xaxes(automargin=True, row=3)
    # and the same for the counts down the left, which the strips have pushed off the figure
    figure.update_yaxes(automargin=True, row=1, col=1)


def add_niche_band(figure, hosts, labelled):
    '''
    The strip of niche under each column -- whether the parasite sits inside a host cell or
    outside it -- and the two legends a figure needs once it carries one, the colours of the
    bars above and the colours of the strip below being two keys to two different parts of
    it.

    Drawn for the figures whose bars are already spending their colour on something else.
    It says which host proteins the parasite is in a position to reach at all, which is a
    fact about the parasite and not about the quantity the bars are drawn from, so it
    belongs under them rather than in them.

    :param figure: a figure host_columns built with a band to spare, modified in place
    :param list hosts: the (host, its rows, its parasites in axis order) of host_columns
    :param set labelled: values already in the legend, added to as they are drawn
    '''
    add_band(figure, hosts, 'niche', web_utils.NICHE_COLORS, labelled,
             unknown=web_utils.NICHE_COLORS[web_utils.UNKNOWN_NICHE])
    figure.update_layout(margin=dict(t=125),
                         legend=dict(y=1.28, title_text='taxonomic group',
                                     title_font=dict(size=11)),
                         legend2=dict(orientation='h', yanchor='bottom', y=1.14, x=0,
                                      title_text=web_utils.NICHE_TITLE,
                                      title_font=dict(size=11), font=dict(size=11)))
    # the names belong under the strip, which is the foot of the figure now, and a column
    # labelled twice is a column labelled once too often. The room they need is taken off
    # the figure rather than out of the margin, which is where the strip has pushed them
    figure.update_xaxes(showticklabels=False, row=1)
    figure.update_xaxes(automargin=True, row=2)


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


def classify_side(df, side):
    '''
    Which localization class DeepLoc puts the protein of `side` in, read on the classes
    that side was filtered on.

    The parasite side went through the secretome filter, which reads the surface pair
    alone, so every parasite protein is read on the same two classes. The host side went
    through apply_deeploc_filter, which reads the classes the niche of the parasite allowed
    -- the surface pair for every parasite, the cytosol and the nucleus for the
    intracellular ones on top -- so a host protein is classified once per parasite reaching
    it and not once for itself: the same membrane protein is a cell membrane protein under
    every parasite, and a cytosolic one as well only under the parasites let into the
    cytosol.

    :param df: rows carrying the DeepLoc probabilities and, for the host side, the `niche`
               of the parasite each row belongs to
    :param str side: 'source' for the parasite protein, 'target' for the host protein
    :return: series of class names, aligned to the rows given
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

    One row per interaction and not per protein, which is what makes this the same quantity
    the bars already draw: the column keeps its height and only gains a split. The protein
    tables below count a protein once however many interactions it is in, and answer a
    different question -- what a parasite reaches -- from this one, which is where the
    predictions themselves land.

    An interaction whose protein the DeepLoc table does not carry keeps its place in the
    count as NOT_SURFACE rather than being dropped, so the bars stay the counts of the
    figure above them whatever the localisations cover.

    :param df: overview predictions, as get_overview_predictions builds them
    :param str data_dir: directory holding the localisations
    :param str side: 'source' for the parasite proteins, 'target' for the host proteins
    :return: the predictions with a `surface` column, or None without a localisation table
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
    about where each protein sits: the probability of each localization class the side was
    filtered on, and which of them the protein was called for -- one class, several of
    them, or none. A host protein is read once per parasite reaching it rather than once
    for itself, since the classes it is read on are the ones that parasite's niche allowed.

    `side` is 'source' for the parasite proteins each parasite reaches its host with, and
    'target' for the host proteins they reach. Both sides went through a localisation
    filter to get here, but not the same one, and so they are not read on the same classes:
    the parasite side through the secretome filter, which allows a multicellular parasite
    nothing but secreted proteins and reads the surface pair alone, and the host side
    through apply_deeploc_filter, which reads all four -- the surface pair for every
    parasite, the cytosol and the nucleus for the intracellular ones on top. Only the host
    side can therefore be read as a comparison between parasites, and there the split of a
    column is partly the niche of the parasite it stands under.

    One row per protein and not per interaction: a parasite protein reaching eleven host
    proteins is one protein, not eleven. A parasite infecting two hosts has its proteins
    counted under each, since the page compares the hosts.

    :param str data_dir: directory holding predictions.parquet and the localisations
    :param dict config: parsed configuration
    :param str side: 'source' for the parasite proteins, 'target' for the host proteins
    :return: one row per host, parasite and protein, or None without a localisation table
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
    The spread of the probability itself, before it is a class: one box per parasite over
    the proteins of one surface class, in the same columns and colours as the figures above
    it. Read within a class rather than across the whole scale, a box says how sure DeepLoc
    was of the probability for the localization it assigned.

    The dotted line is the cut-off of the class, which is what the secretome filter kept
    these proteins on, so nothing sits below it and the scale starts just under it: a box
    resting on the line is a parasite whose proteins only just qualified, a box near 1 one
    whose proteins are unambiguous.

    A parasite with no protein of the class has no box. The membrane figure is handed the
    unicellular parasites alone: the secretome filter keeps a multicellular parasite's
    proteins for being secreted and for nothing else, so the handful of them DeepLoc also
    assigns to the cell membrane are an accident of which secreted proteins carry a second
    assignment rather than a sample of anything, and a box over one such protein says
    nothing about the parasite it is drawn under.

    Every protein is drawn as well as summarised, jittered across the width of its box.
    Several of these boxes stand on ten or twenty proteins, which is too few for a box to
    be read as a distribution without seeing what is behind it. The confidence figure at
    the top of the page keeps its points off for the opposite reason: there a point is an
    interaction rather than a protein, and there are thousands of them to a column.

    :param proteins: the proteins of one surface class, as get_parasite_proteins builds them
    :param dict palette: {taxonomic group: colour}
    :param float width: pixels the figure is drawn across, for the names over the columns
    :param str score: the probability column to draw
    :param float cutoff: the cut-off of that class, drawn as a dotted line
    :param str y_title: what the probability is called down the left of the figure
    :param float point_size: diameter of a protein drawn beside its box, in pixels. The
                             proteins are jittered across the width of the box and the box
                             itself is left unfilled so they are read through it; a column
                             narrower than its proteins are many comes out as a cloud rather
                             than as points that can be counted
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
    # the cut-off belongs to the boxes and not to the strip under them, so it is drawn in
    # the row of the columns rather than across the whole figure
    figure.add_hline(y=cutoff, line_width=1, line_dash='dot', line_color='#969696',
                     row=1, col='all')
    add_niche_band(figure, hosts, set())

    return figure


@st.cache_data(show_spinner=False)
def generate_host_score_boxes(proteins, width, point_size=3):
    '''
    How sure DeepLoc was of the host proteins, one column per localization class and one
    box per host inside it.

    Grouped by class and not by host: the four classes are called on four different
    probabilities and at four different cut-offs, so what a box means changes with the class
    and not with the host. A column per class puts the hosts on a common axis and leaves one
    cut-off to draw per column instead of four in every column.

    Every host protein is above the cut-off of the class it is drawn under, the host filter
    reading those same thresholds, so the scale starts just under the lowest of them.

    Four columns and not one per combination of classes. A protein DeepLoc calls for more
    than one class is drawn in each of them, at that class's own probability: the
    probabilities are separate statements about the protein and each belongs under the class
    it is about, where a column of its own would have to pick one of them to stand for the
    rest. It is also how the parasite figures below read their proteins.

    The cytosol and nucleus columns stand on the host proteins of the intracellular
    parasites alone -- no other parasite was allowed one -- so they are a statement about
    those hosts under those parasites and not about the host pool as a whole.

    Per host and not per parasite. The parasites of a host draw their interactors from the
    same few hundred host proteins, so a box per parasite is a box over a sample of one
    pool: their medians span 0.65 to 0.79 with no order to them, forty-six near-copies of
    the pool and of each other.

    A host protein is counted once however many parasites reach it, since here it is a
    protein of the host and not an interactor of anything.

    :param proteins: the host proteins, as get_interactor_proteins(side='target') builds them
    :param float width: pixels the figure is drawn across, for the names over the columns
    :param float point_size: diameter of a protein drawn beside its box, in pixels
    '''
    # one row per protein and class it was called for, so a protein over two cut-offs is a
    # row under each. Read from the probability rather than from the class of the protein,
    # which names one class for a protein called for several; proteins over no cut-off at
    # all are left out by there being no column they belong under. A class the table
    # carries no probability for is no column either, which is what a snapshot data
    # directory written before those columns existed leaves out
    columns = []
    for surface_class in web_utils.HOST_CLASSES:
        score_column = web_utils.DEEPLOC_SCORES[surface_class]
        if score_column not in proteins:
            continue
        # only the proteins some parasite was allowed to meet in that class: a membrane
        # protein of an extracellular parasite is often cytosolic too, and drawing it in
        # the cytosol column would put a protein there that nothing reaches there
        allowed = proteins[proteins['niche'].map(
            lambda niche: surface_class in web_utils.niche_classes(niche))]
        counted = allowed.drop_duplicates(['host', 'protein'])
        in_class = counted[counted[score_column] > web_utils.DEEPLOC_CUTOFFS[surface_class]]
        columns.append(in_class.assign(surface=surface_class, score=in_class[score_column]))
    scored = pd.concat(columns, ignore_index=True)
    # host_columns splits on 'host' and orders what is in a column by 'group_rank'; here a
    # column is a localization class and what is in it are the hosts, so the two are
    # swapped. The columns come out in the order the frames were concatenated, which is the
    # order HOST_CLASSES declares and the order the figures above split their bars in
    host_order = {host: rank for rank, host in enumerate(scored['host'].unique())}
    # the species abbreviated the way the columns of the figures above abbreviate a
    # parasite: four classes to a figure leave a quarter of the room the two of them left,
    # and four full names to a column are written over each other
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
        # the cut-off this class was called at, drawn in the column it applies to rather
        # than across the figure, where it would read as a threshold on both
        figure.add_hline(y=web_utils.DEEPLOC_CUTOFFS[surface_class], line_width=1,
                         line_dash='dot', line_color=SURFACE_LINE_COLORS[surface_class],
                         row=1, col=column)

    figure = style_host_columns(figure, 'P(the assigned localization)')
    # the scale clears the cut-offs of the classes actually drawn, not of every class there
    # is, so a directory carrying two of them is not given the headroom of four
    drawn = [web_utils.DEEPLOC_CUTOFFS[c] for c in scored['surface'].unique()]
    figure.update_yaxes(range=[score_floor(*drawn, scored['score'].min()), 1],
                        automargin=True, row=1, col=1)
    figure.update_xaxes(tickangle=0, automargin=True)

    return figure


@st.cache_data(show_spinner=False)
def generate_interactions_per_parasite(df, palette, width, score, classes=None, share=False):
    '''
    How many interactions are predicted for each parasite at or above a confidence, in one
    column per host, split by where DeepLoc puts the protein on one side of each of them.

    The columns are laid out from every prediction and only the bars are thresholded, so a
    parasite left with nothing keeps its place on the axis and is read as an absence rather
    than disappearing out of a figure the one below it still draws it in.

    A bar is still the interactions of a parasite and the stack is the same interactions
    divided up, so the figure is read two ways at once: how much a parasite is predicted to
    interact, and where those interactions land. `share` normalises each column to itself
    and asks the second question alone -- the counts run from 39 interactions to two and a
    half thousand, and a segment of the smallest column is invisible beside the largest.

    On the host side the split is read against the niche in the lower strip, and the two say
    different things. The niche is what the filter allowed: the oranges can only appear
    under an intracellular parasite, apply_deeploc_filter keeping the cytosol and the
    nucleus for those alone. How much orange there is, and how the blue divides, is then a
    difference between parasites that were allowed the same classes.

    :param df: overview predictions, with a `surface` column where the bars are to be split
    :param dict palette: {taxonomic group: colour}
    :param float width: pixels the figure is drawn across, for the names over the columns
    :param float score: the confidence to count from
    :param classes: the classes to stack, in the order they are stacked in, or None to
                    colour whole bars by taxonomic group instead -- which is what a data
                    directory without a localisation table is left with
    :param bool share: draw each column as its own 100% rather than as a count
    '''
    if classes is None:
        return generate_interactions_by_group(df, palette, width, score)

    # a class nothing was called for is still an entry of the legend, so the two sides are
    # read on the classes their filter allowed and not on the ones that happen to be full.
    # NOT_SURFACE is not one of them: it is a protein the localisations do not carry, kept
    # in the count so the column stays the count, and drawn only where there is one
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
        # the whole column, which a segment is a share of and which the hover names: a
        # parasite the threshold emptied has none, and is drawn as the absence it is
        total = counts.sum(axis=1).replace(0, pd.NA)
        for surface_class in classes:
            figure.add_trace(
                go.Bar(x=counts.index, y=counts[surface_class] / total if share
                                         else counts[surface_class],
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

    figure = style_host_columns(figure, 'share of predicted interactions' if share
                                        else 'predicted interactions')
    figure.update_layout(barmode='stack', bargap=0.2)
    if share:
        figure.update_yaxes(range=[0, 1], tickformat='.0%', row=1)
    add_group_and_niche_bands(figure, hosts, palette)

    return figure


def generate_interactions_by_group(df, palette, width, score):
    '''
    The same counts with nothing to split them by: whole bars in the colour of the parasite's
    taxonomic group, with the niche in the one strip under them. What a data directory
    without a DeepLoc table is drawn as.

    :param df: overview predictions, as get_overview_predictions builds them
    :param dict palette: {taxonomic group: colour}
    :param float width: pixels the figure is drawn across, for the names over the columns
    :param float score: the confidence to count from
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
                       # a bar per parasite, whatever group it is drawn in the colour of,
                       # so the bars stand over the segments of the strip below them
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

    The one figure of the page the slider does not filter, drawn on every prediction with
    the slider as a line across it instead. Filtering a distribution by the axis it is
    drawn on says nothing -- every lower whisker becomes the slider -- where the line makes
    this figure the key to the slider: how much of a parasite's box stands above it is how
    much of that parasite survives into the counts above.

    The outlying points are left off. Every interaction is a point, and a few hundred of
    them beside each box hide the boxes themselves.

    :param df: overview predictions, as get_overview_predictions builds them
    :param dict palette: {taxonomic group: colour}
    :param float width: pixels the figure is drawn across, for the names over the columns
    :param float score: where the counts above are cut, drawn as a line
    '''
    figure, hosts = host_columns(df, width, bands=1)
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
    # the line belongs to the boxes and not to the strip under them, so it is drawn in the
    # row of the columns rather than across the whole figure
    figure.add_hline(y=score, line_width=1, line_dash='dot', line_color='#969696',
                     row=1, col='all')
    # the line is drawn in every column and named in the first: an annotation per column is
    # the same three words written once for each host
    figure.add_annotation(text='counted above', x=0, y=score, xref='x domain', yref='y',
                          xanchor='left', yanchor='bottom', showarrow=False,
                          font=dict(size=10, color='#969696'), row=1, col=1)
    add_niche_band(figure, hosts, set())

    return figure


st.caption('Protein-protein interactions between parasites and their hosts, predicted by '
           'orthology transfer and restricted to the host proteins expressed in a tissue '
           'the parasite is known to infect. This page presents every prediction per host; '
           '**Parasites of a host** compares the parasites of a single host, **Hosts of a '
           'parasite** follows one parasite across the hosts it infects, and **Host-parasite '
           'network** shows the network of one host-parasite pair.')
st.markdown("---")

# the overview of what is in the database, first since it is what every other figure is
# drawn on
st.subheader('Hosts and parasites')
hosts_parasites.show(config)
st.markdown("---")

# the figures are stretched to the page, and the names over their columns have to be
# written to fit the column each of them came out with, so the page is measured once here
# and every figure drawn to the width it reports. Nothing waits on it: page_width answers
# with a laptop until the browser has replied, and the figures are redrawn on the run it does
page = web_utils.page_width()
overview = get_overview_predictions(data_dir, config)
parasite_palette = config.get('parasite_groups', {})
coverage = host_coverage_caption(data_dir, config)
if coverage:
    st.caption(coverage)

# one slider for the counts and the proportions of the page, so what is counted here is
# what the network page opens on rather than several times more of it. It sits in the
# middle of three columns, as the sliders of the other pages do: left to itself a slider
# takes the whole width of the page, which is a metre of track for a range of half a point
with st.columns(3)[1]:
    score = st.slider('Confidence score', MIN_SCORE, MAX_SCORE, DEFAULT_SCORE,
                      help='Interactions predicted below this confidence are left out '
                           'of the counts and of the localization split drawn on them. The '
                           'boxplots keep every prediction: the confidence figure draws '
                           'this threshold as a line instead, and the localization ones '
                           'stand on too few proteins to be thresholded as well.')

st.subheader("Number of predicted interactions per parasite")
# the two toggles over the figure: which end of the interaction the bars are split by, and
# whether the columns are counts or shares. Narrow columns, since a segmented control left
# to itself is stretched over the width of the page
controls = st.columns([1.4, 1, 1.6])
with controls[0]:
    split_side = st.segmented_control(
        'Colour by the localization of', list(SPLIT_SIDES), default=list(SPLIT_SIDES)[0],
        key='split_side',
        help='Every prediction is between one parasite protein and one host protein, and '
             'the bars can be split by where DeepLoc 2 puts either of them.')
with controls[1]:
    bar_scale = st.segmented_control(
        'Bars show', BAR_SCALES, default=BAR_SCALES[0], key='bar_scale',
        help='Counts are the predicted interactions themselves; shares divide each column '
             'by its own total, which is how the split of a parasite with few interactions '
             'is compared with one that has thousands.')
side, split_classes = SPLIT_SIDES[split_side or list(SPLIT_SIDES)[0]]
interactions = get_interaction_localisations(overview, data_dir, side)

if interactions is None:
    st.caption('Predicted interactions per parasite at or above the confidence set above, '
               'grouped by host and coloured by parasite taxonomic group. ' + NICHE_STRIP)
elif side == 'target':
    st.caption('Predicted interactions per parasite at or above the confidence set above, '
               'grouped by host and split by the subcellular localization DeepLoc 2 '
               'predicts for the host protein of each interaction, over the four classes '
               'the host filter reads — extracellular, cell membrane, cytoplasm and '
               'nucleus. Each parasite is read on the classes its niche let the filter keep '
               'a host protein for: the surface pair for every parasite, the cytosol and '
               'the nucleus for the ones with an intracellular stage, which is why the two '
               'oranges appear under those alone. A host protein called for more than one '
               'of the classes its parasite can reach is counted as several. ' + BANDS_STRIP)
else:
    st.caption('Predicted interactions per parasite at or above the confidence set above, '
               'grouped by host and split by the localization DeepLoc 2 assigns the '
               'parasite protein of each interaction — cell membrane, extracellular, '
               'or both. The secretome filter admits a multicellular parasite nothing but '
               'its secreted proteins, so under the cestodes, nematodes and trematodes the '
               'split is the filter rather than the parasite; only the unicellular '
               'parasites had both classes open to them. ' + BANDS_STRIP)

st.plotly_chart(
    generate_interactions_per_parasite(overview if interactions is None else interactions,
                                       parasite_palette, page, score,
                                       classes=None if interactions is None else split_classes,
                                       share=bar_scale == BAR_SCALES[1]),
    width='stretch')

st.subheader("Interaction confidence scores")
st.caption('Boxplots of the distribution of confidence scores per parasite. Each score is '
           'the average of the experimental and database evidence channels of the KOG-KOG '
           'link the prediction was transferred from. Every prediction is counted here, '
           'whatever the slider is set to; the dotted line marks the threshold it selects '
           'instead. ' + NICHE_STRIP)
st.plotly_chart(generate_confidence_per_parasite(overview, parasite_palette, page, score),
                width='stretch')

# the boxes of the probability itself, one protein per row rather than one interaction:
# they are asking how sure DeepLoc was of a protein, which the interactions it turns up in
# do not bear on. Drawn on every prediction and not at the threshold, since at 0.7 a
# parasite is left a median of eight of its own proteins, which is not a distribution
host_proteins = get_interactor_proteins(data_dir, config, 'target')
parasite_proteins = get_interactor_proteins(data_dir, config, 'source')
# the membrane figure is drawn over the unicellular parasites alone
unicellular = None
if parasite_proteins is not None:
    unicellular = parasite_proteins[parasite_proteins['group'].isin(UNICELLULAR_GROUPS)]

if host_proteins is not None:
    st.subheader("Localization probabilities of host proteins")
    st.caption('Boxplots of the DeepLoc 2 probabilities of the host proteins for their '
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
    st.plotly_chart(generate_host_score_boxes(host_proteins, page), width='stretch')

if parasite_proteins is not None:
    st.subheader("Localization probabilities of extracellular parasite proteins")
    st.caption('Boxplots of the DeepLoc 2 extracellular probability for proteins assigned '
               'extracellular or both classes for each parasite, over every prediction whatever '
               'the slider is set to. The dotted line marks the cut-off, DeepLoc 2\'s default '
               'threshold for extracellular under the Accurate model '
               f'({web_utils.DEEPLOC_CUTOFFS[web_utils.EXTRACELLULAR]:.3f}), which is '
               'what the secretome filter kept these proteins on. Individual proteins are '
               'shown as points behind each box; for parasites with a hundred or more proteins '
               'the points are read as density. ' + NICHE_STRIP)
    st.plotly_chart(
        generate_surface_scores_per_parasite(
            parasite_proteins[parasite_proteins['surface'].isin(
                [web_utils.EXTRACELLULAR, web_utils.BOTH_SURFACE])],
            parasite_palette, page, 'extracellular',
            web_utils.DEEPLOC_CUTOFFS[web_utils.EXTRACELLULAR],
            # forty-five columns to a row and up to a hundred and eighty proteins in one of
            # them, so the smallest dot that still carries colour
            'P(extracellular)', point_size=2),
        width='stretch')

    st.subheader("Localization probabilities of membrane parasite proteins")
    st.caption('The equivalent for parasite proteins assigned to cell membrane or both '
               'classes, scored on that probability, again over every prediction. Only '
               'unicellular parasites are represented, the secretome filter admitting a '
               'multicellular parasite nothing but its secreted proteins. Individual '
               'proteins are shown as points; boxes over very few proteins (one each for '
               '*G. lamblia* and *V. corneae*) should not be read as distributions. The '
               'dotted line marks the cut-off, DeepLoc 2\'s default for cell membrane under '
               f'the Accurate model ({web_utils.DEEPLOC_CUTOFFS[web_utils.CELL_MEMBRANE]:.3f}). '
               + NICHE_STRIP)
    st.plotly_chart(
        generate_surface_scores_per_parasite(
            unicellular[unicellular['surface'].isin(
                [web_utils.CELL_MEMBRANE, web_utils.BOTH_SURFACE])],
            parasite_palette, page, 'cell_membrane',
            web_utils.DEEPLOC_CUTOFFS[web_utils.CELL_MEMBRANE],
            # a third of the parasites and a fifth of the proteins of the figure above, so
            # the columns are wide enough for the proteins to be told apart
            'P(cell membrane)', point_size=3),
        width='stretch')

st.markdown("---")


# Footer
with st.container():
    web_utils.footer()
