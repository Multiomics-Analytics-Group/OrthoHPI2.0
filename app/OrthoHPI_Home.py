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

config = utils.read_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()

# fallback for a parasite without a `group` in the config
UNKNOWN_GROUP = 'Unclassified'
# the membrane figure is drawn over these: a multicellular parasite reaches its host with
# secreted proteins alone
UNICELLULAR_GROUPS = ('Apicomplexa', 'Kinetoplastida', 'Other protozoa', 'Microsporidia')
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
# same default and range as the network page, so the counts agree with the network opened
# next
MIN_SCORE, MAX_SCORE, DEFAULT_SCORE = 0.35, 0.9, 0.35
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
# which end of an interaction the counts figure is split by; host side first
SPLIT_SIDES = {'Host proteins': ('target', HOST_SPLIT_CLASSES),
               'Parasite proteins': ('source', PARASITE_SPLIT_CLASSES)}
# the counts run from 39 to two and a half thousand, so the split of the small parasites is
# only legible in shares
BAR_SCALES = ('Counts', 'Share')
# the whole proteome is drawn this faint behind the solid bar of its eligible subset
PROTEOME_OPACITY = 0.3
# how far under the lowest cut-off or point a probability scale starts
SCALE_MARGIN = 0.05
# rows of the shared-family dot plot, pixels a row takes, and what the strip, the names
# and the legends take besides
TOP_FAMILIES = 40
DOT_ROW = 19
DOT_CHROME = 230
# most gene symbols a family is named by before the rest are left to the hover, and the
# most characters
SYMBOLS_IN_LABEL = 3
LABEL_CHARS = 24
# room a family name takes left of the dots: ~6.2 px a character, plus the axis title
LABEL_CHAR = 6.2
LABEL_PADDING = 40
# height of one horizontal legend row, in pixels
LEGEND_ROW = 46


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
    # the order of the circos: taxonomic group as configured, then name, unclassified last
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
        parasites = host_df[['taxid1_label', 'name', 'group_rank']].drop_duplicates()
        parasites = parasites.sort_values(by=['group_rank', 'taxid1_label'], kind='stable')
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
        segments = dict(tuple(host_df.groupby(field, observed=True)))
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
def generate_interactions_per_parasite(df, palette, width, score, classes=None, share=False):
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
def generate_proteome_sizes(config, sizes, eligible, palette):
    '''
    One bar per parasite: the whole proteome STRING holds for it, faint, with the
    proteins the secretome filter let through solid in front. The ratio is what the
    counts above stand on -- a small proteome, or a multicellular parasite kept to its
    secreted proteins, has few interactions to offer before any prediction is made.
    '''
    order = {g: i for i, g in enumerate(palette)}
    pool = eligible['taxid'].astype(str).value_counts()
    sizes = sizes.set_index('taxid')['proteins']
    rows = []
    for taxid, parasite in config['parasites'].items():
        group = parasite.get('group', UNKNOWN_GROUP)
        rows.append({'name': short_name(parasite['label']), 'group': group,
                     'group_rank': order.get(group, len(order)),
                     'label': parasite['label'],
                     # what the secretome filter admits of each: a multicellular parasite
                     # reaches its host with secreted proteins alone
                     'kept': 'secreted' if parasite.get('multicellular')
                             else 'secreted or membrane',
                     'proteome': sizes.get(str(taxid)), 'eligible': pool.get(str(taxid), 0)})
    df = pd.DataFrame(rows).dropna(subset=['proteome'])
    df = df.sort_values(by=['group_rank', 'label'], kind='stable')
    df['share'] = df['eligible'] / df['proteome']

    figure = go.Figure()
    for group, group_df in df.groupby('group', sort=False):
        color = palette.get(group, UNKNOWN_COLOR)
        custom = group_df[['proteome', 'eligible', 'share', 'kept']]
        hover = ('%{x}<br>%{customdata[0]:,} proteins in STRING<br>%{customdata[1]:,} '
                 '%{customdata[3]} (%{customdata[2]:.0%})' f'<extra>{group}</extra>')
        figure.add_trace(
            go.Bar(x=group_df['name'], y=group_df['proteome'], name=group,
                   marker_color=color, opacity=PROTEOME_OPACITY, customdata=custom,
                   legendgroup=group, showlegend=False, hovertemplate=hover))
        figure.add_trace(
            go.Bar(x=group_df['name'], y=group_df['eligible'], name=group,
                   marker_color=color, customdata=custom,
                   legendgroup=group, hovertemplate=hover))

    figure.update_layout(barmode='overlay', bargap=0.2, height=420, plot_bgcolor='white',
                         margin=dict(l=0, r=0, t=40, b=10),
                         legend=dict(orientation='h', yanchor='bottom', y=1.02, x=0,
                                     title_text='', font=dict(size=11)))
    figure.update_xaxes(categoryorder='array', categoryarray=df['name'].tolist(),
                        tickangle=-60, showgrid=False, tickfont=dict(size=11),
                        automargin=True)
    figure.update_yaxes(title_text='proteins', showgrid=True, gridcolor='#f0f0f0',
                        zerolinecolor='#e0e0e0', automargin=True)

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
    # in the row of the columns rather than across the strip under them
    figure.add_hline(y=score, line_width=1, line_dash='dot', line_color='#969696',
                     row=1, col='all')
    # the line is drawn in every column and named in the first
    figure.add_annotation(text='counted above', x=0, y=score, xref='x domain', yref='y',
                          xanchor='left', yanchor='bottom', showarrow=False,
                          font=dict(size=10, color='#969696'), row=1, col=1)
    add_niche_band(figure, hosts, set())

    return figure


def family_label(names):
    '''
    The name a host family is drawn under: the gene symbols of its proteins, since the
    orthology group id names nothing to read. The symbols are upper-cased so the same
    family reads the same way whichever hosts it is reached in.
    '''
    names = sorted(set(str(n).upper() for n in names))
    named = []
    for name in names[:SYMBOLS_IN_LABEL]:
        # the first name goes in whatever its length, so a family is never drawn under an
        # ellipsis alone
        if named and len(', '.join(named + [name])) > LABEL_CHARS:
            break
        named.append(name)
    label = ', '.join(named)

    return f'{label}…' if len(named) < len(names) else label


@st.cache_data(show_spinner=False)
def get_top_shared_families(df, score, top=TOP_FAMILIES):
    '''
    The host protein families the most parasites are predicted to interact with, across
    every host: a family is the orthology group of the host protein, which is what the
    same protein of two hosts has in common. One row per parasite, host and family.
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

    families = list(ranked.head(top).index)
    dots = dots[dots['group2'].isin(families)].copy()
    dots['parasites'] = dots['group2'].map(ranked['parasites'])
    dots['pairs'] = dots['group2'].map(ranked['pairs'])
    dots['host proteins'] = dots['targets'].map(', '.join)
    # named by every symbol the family carries in any host, not only the ones this dot
    # reaches
    labels = {family: family_label(n for names in rows['targets'] for n in names)
              for family, rows in dots.groupby('group2')}
    dots['family'] = dots['group2'].map(labels)
    # the columns of the skeleton: the same parasite annotations the other figures carry
    parasite_rows = kept[['host', 'taxid1_label', 'name', 'group', 'group_rank',
                          'niche']].drop_duplicates()
    dots = dots.merge(parasite_rows, on=['host', 'taxid1_label'])

    return dots, [labels[f] for f in families], len(ranked)


@st.cache_data(show_spinner=False)
def generate_shared_family_dots(df, dots, families, palette, width):
    '''
    A dot wherever a parasite is predicted to interact with a protein of one of the
    families, in the columns the figures above use, the families ordered by how many
    parasites reach them so the most widely reached is the top row.
    '''
    height = DOT_ROW * len(families) + DOT_CHROME
    # the strip keeps the height it has under the other figures
    figure, hosts = host_columns(df, width, bands=1,
                                 band_height=BAND_HEIGHT * 470 / height)
    rows_of = {f: i for i, f in enumerate(families)}
    dots = dots.assign(y=dots['family'].map(rows_of))
    # the largest dot is as wide as a column has room for, its area standing for the
    # largest degree
    columns = sum(max(len(names), MIN_COLUMN) for _, _, names in hosts)
    room = width - (LABEL_CHAR * max(len(f) for f in families) + LABEL_PADDING)
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
                           customdata=rows[['family', 'parasites', 'pairs', 'degree',
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

    figure = style_host_columns(figure, 'host protein family')
    add_niche_band(figure, hosts, labelled)
    # reversed, so the most-shared family is the top row; the ticks are set on every
    # column so the grids line up
    figure.update_yaxes(range=[len(families) - 0.5, -0.5], tickmode='array',
                        tickvals=list(range(len(families))), ticktext=families,
                        zeroline=False, row=1)
    # the names are longer than the counts the margin was set for
    figure.update_yaxes(automargin=True, row=1, col=1)
    figure.update_xaxes(showgrid=True, gridcolor='#f0f0f0', row=1)
    # the legends are placed from the top of the figure rather than the plot, which the
    # rows have made tall
    figure.update_layout(height=height,
                         legend=dict(yref='container', yanchor='top', y=1),
                         legend2=dict(yref='container', yanchor='top',
                                      y=1 - LEGEND_ROW / height))

    return figure


st.caption('Protein-protein interactions between parasites and their hosts, predicted by '
           'orthology transfer and restricted to the host proteins expressed in a tissue '
           'the parasite is known to infect. This page presents every prediction per host; '
           '**Parasites of a host** compares the parasites of a single host, **Hosts of a '
           'parasite** follows one parasite across the hosts it infects, and **Host-parasite '
           'network** shows the network of one host-parasite pair.')
st.markdown("---")

st.subheader('Hosts and parasites')
overview = get_overview_predictions(data_dir, config)
hosts_parasites.show(config, interactions=len(overview))
st.markdown("---")

# measured once here so the column titles can be fitted; page_width answers with a default
# until the browser has replied
page = web_utils.page_width()
parasite_palette = config.get('parasite_groups', {})
coverage = host_coverage_caption(data_dir, config)
if coverage:
    st.caption(coverage)

# one slider for the page, in the middle column so it does not take the whole width
with st.columns(3)[1]:
    score = st.slider('Confidence score', MIN_SCORE, MAX_SCORE, DEFAULT_SCORE,
                      help='Interactions predicted below this confidence are left out '
                           'of the counts and of the localization split drawn on them. The '
                           'boxplots keep every prediction: the confidence figure draws '
                           'this threshold as a line instead, and the localization ones '
                           'stand on too few proteins to be thresholded as well.')
    # what the threshold keeps, under the slider
    kept = int((overview['weight'] >= score).sum())
    st.caption(f'{kept:,} of {len(overview):,} interactions ({kept / len(overview):.0%}) '
               f'at or above {score:g}.')

st.subheader("Number of predicted interactions per parasite")
# narrow columns: a segmented control left to itself is stretched over the page
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

proteome_sizes = web_utils.load_proteome_sizes(data_dir)
eligible_proteins = web_utils.load_eligible_proteins(data_dir)
if proteome_sizes is not None and eligible_proteins is not None:
    with st.expander('Parasite proteome sizes'):
        st.caption('The whole proteome of each parasite as STRING holds it, faint, and in '
                   'front of it the proteins the secretome filter let through — the '
                   'membrane and secreted proteins of a unicellular parasite, the secreted '
                   'proteins alone of a multicellular one.')
        st.plotly_chart(generate_proteome_sizes(config, proteome_sizes, eligible_proteins,
                                                parasite_palette),
                        width='stretch')

st.subheader("Host protein families common to several parasites")
shared_families = get_top_shared_families(overview, score)
if shared_families is None:
    st.info('No host protein family is reached by more than one parasite at this confidence.')
else:
    family_dots, families, shareable = shared_families
    shown = len(families)
    selection = (f'The {shown} host protein families reached by the most parasites, of the '
                 f'{shareable} reached by more than one, ' if shareable > shown else
                 f'The {shown} host protein families reached by more than one parasite, ')
    st.caption(selection + 'at or above the confidence set above. A family is the '
               'orthology group of the host protein, and is named by the gene symbols of its proteins. A dot '
               'wherever a parasite is predicted to interact with a protein of the family, '
               "sized by the number of that parasite's proteins reaching it; hover a dot for the host proteins behind it. "
               + NICHE_STRIP)
    st.plotly_chart(generate_shared_family_dots(overview, family_dots, families,
                                                parasite_palette, page),
                    width='stretch')

st.subheader("Interaction confidence scores")
st.caption('Boxplots of the distribution of confidence scores per parasite. Each score is '
           'the average of the experimental and database evidence channels of the KOG-KOG '
           'link the prediction was transferred from. Every prediction is counted here, '
           'whatever the slider is set to; the dotted line marks the threshold it selects '
           'instead. ' + NICHE_STRIP)
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
            # forty-five columns and up to 180 proteins in one, so the smallest dot that
            # still carries colour
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
            # fewer parasites and proteins than the figure above, so the dots can be told
            # apart
            'P(cell membrane)', point_size=3),
        width='stretch')

st.markdown("---")


# Footer
with st.container():
    web_utils.footer()
