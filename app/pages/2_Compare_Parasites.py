import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import utils
import web_utils
import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from css import style

st.set_page_config(layout="wide", page_title="OrthoHPI 2.0", menu_items={})
style.load_css()
web_utils.show_header('Parasites of a host')

# Read dataset
config = utils.read_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()


# fallback for a parasite without a `group` in the config
UNKNOWN_GROUP = 'Unclassified'
UNKNOWN_COLOR = '#999999'
# confidence the figures of shared interactors start at, and the range the slider spans.
# The same default and range as the network page, so a parasite shows the same
# interactions here as it does there
MIN_SCORE, MAX_SCORE, DEFAULT_SCORE = 0.35, 0.9, 0.35
# side of a cell of the shared-interactor heatmap, in pixels, which is what sizes that figure
CELL = 22
# the diagonal of that heatmap, where a parasite meets itself: a mid grey, so that it reads
# as a boundary between the two halves rather than as a count on the blue scale
DIAGONAL_COLOUR = '#9e9e9e'
# the room the colour bar of that heatmap takes to the right of the matrix, in pixels: the
# gap it is set out from the matrix by, then the bar, the counts written beside it -- about
# six and a half pixels a digit at ten points -- and its rotated title. The gap is held in
# pixels, rather than the fraction of the plot a colour bar is set out by by default, so
# that the room the bar needs is the same on every screen: plotly widens a margin that is
# short of what sits in it, and a widened right margin is a plot area no longer square
COLORBAR_GAP = 8
COLORBAR_ROOM = 71
COLORBAR_DIGIT = 6.5
# how a similarity is written wherever the figure gives one: beside the colour bar, where
# it sets the room the bar needs, and in the hover of a cell. Two decimals: most pairs of a
# host sit under a tenth -- the median pair of a human is 0.08 -- and written to one they
# are all the same number
TICK = '.2f'
# the height of one row of a legend below that heatmap, in pixels, which is what the room
# left under the matrix is counted in
LEGEND_ROW = 22
# the room one entry of a legend takes across the foot of that heatmap, in pixels:
# its marker and the padding around it, and about seven and a half pixels a character of the
# longest name in it -- plotly gives every entry the room the longest of them needs
LEGEND_ENTRY = 40
LEGEND_CHAR = 7.5
# the room left under the last legend of a dot matrix, in pixels
LEGEND_PAD = 14
# and about six pixels a character of the longest name in a legend of one, which is the
# room plotly gives every entry of it. Narrower than the LEGEND_CHAR of the heatmap: the
# entries there are set out to the width of the longest, the eleven-point names here pack
# closer than that, and a legend measured too wide is a blank row between two legends
DOT_LEGEND_CHAR = 6.0
# and the room the title of one takes above its entries. A horizontal legend is titled on
# top -- plotly defaults `title.side` to that on an orientation of `h` -- so the title is a
# row of the legend and not a first entry set beside the others
LEGEND_TITLE_ROW = 18
# the smallest the parasite names on the axis of a matrix are written, in points
SMALLEST_LABEL = 9
# session key the shared-interactor matrix counts its remounts under, which is what lets
# the same cell be opened twice running
CELL_NONCE_KEY = 'shared_cell_nonce'


# The localization of a host protein is a band beside its row of the shared-interactors
# dot plot, in the colours web_utils keys every localization figure of the app with, and
# in this order. It is a property of the protein and so constant along a row, which is
# what makes it a band: as the shape of the dot it was drawn once per parasite reaching
# the protein, at a size where a circle and a diamond cannot be told apart anyway, and it
# distorted the degree the same marker carries in its area.
# NO_LOCALISATION is a protein DeepLoc was never run on, or one whose data directory
# predates pipeline/build_deeploc_localisations.py
NO_LOCALISATION = 'Not available'
LOCALISATION_ORDER = list(web_utils.HOST_CLASSES) + [web_utils.SEVERAL,
                                                     web_utils.NOT_SURFACE, NO_LOCALISATION]
# paler than the grey of a protein the filter had a call for and rejected, a row DeepLoc
# says nothing at all about being the emptier of the two
LOCALISATION_COLORS = {**web_utils.LOCALISATION_COLORS, NO_LOCALISATION: '#f0f0f0'}
# names the DeepLoc columns are read under in the hover of the shared-interactors matrix
DEEPLOC_LABELS = {'surface': 'DeepLoc', 'localizations': 'localizations',
                  **{column: f'P({name.lower()})'
                     for name, column in web_utils.DEEPLOC_SCORES.items()}}


@st.cache_data(show_spinner=False)
def count_interactions_per_tissue(data_dir, config, host_taxids, score=MIN_SCORE):
    '''
    Predicted interactions per parasite and tissue, and per parasite, tissue and cell type,
    keeping only the tissues each parasite is known to infect (config['parasites']).

    `score` drops the interactions predicted below that confidence, the same cut
    get_tissue_expressed_predictions applies, so the figures of this page are all counting
    the same interactions whatever the slider is set to.

    Both are counted as distinct (parasite protein, host protein) pairs at their own level,
    which is the only honest way to size a tissue: a host protein is expressed in several
    cell types of a tissue, so a tissue counted as the sum of its cell types counts the
    same interaction once per cell type. That number is how finely the HPA annotates the
    tissue -- lung has 13 cell types, blood has one -- and not how many interactions can
    take place there. The two frames therefore do not add up to each other, on purpose.

    The tissue counts take a host protein to be present wherever the HPA annotates it,
    while the cell type counts keep only cell types where it exceeds 1 nTPM
    (web_utils.keep_expressed_cell_types). Host proteins the HPA gives no cell type --
    every host but human, the single cell data being human only -- are counted in their
    tissue and left out of the cell type frame.
    '''
    predictions = web_utils.get_host_predictions(data_dir, host_taxids)[
        ['taxid1', 'taxid1_label', 'source', 'target', 'target_name', 'weight']]
    predictions = predictions[predictions['weight'] >= score].drop('weight', axis=1)
    tissues = utils.read_parquet_file(input_file=f'{data_dir}/tissues_cell_types.parquet')
    tissues = tissues.rename({'Gene': 'target'}, axis=1)[['target', 'Tissue', 'Cell type',
                                                         'nTPM']]

    mapped_tissues = config['tissues']
    infected_tissues = pd.DataFrame([(str(taxid), mapped_tissues[t].lower())
                                     for taxid, parasite in config['parasites'].items()
                                     for t in parasite['tissues']],
                                    columns=['taxid1', 'Tissue'])

    aux = predictions.drop_duplicates().astype({'taxid1': str})
    aux = pd.merge(aux, tissues, on='target')
    aux = pd.merge(aux, infected_tissues, on=['taxid1', 'Tissue'])

    # The pair is keyed by the display gene name so aliases of the same host gene do not
    # inflate tissue or cell-type interaction counts.
    def pairs(frame, *level):
        return (frame.drop_duplicates(list(level) + ['source', 'target_name'])
                     .groupby(list(level), observed=True).size()
                     .rename('interactions').reset_index())

    return (pairs(aux, 'taxid1_label', 'Tissue'),
            pairs(web_utils.keep_expressed_cell_types(aux), 'taxid1_label', 'Tissue',
                  'Cell type'))


def dot_size(parasites, labels, column):
    '''
    Diameter of the largest dot of a matrix of `parasites` columns whose rows are named by
    `labels`, drawn in a column `column` pixels wide.

    A column of the matrix is what the figure has left once the row labels have taken the
    room they need, divided between the parasites: forty parasites beside a protein name is
    a column of a few pixels, and a dot drawn at plotly's default fifteen there is a row of
    dots run together into a bar.
    '''
    # 6.2 pixels a character is the width of the default axis font, and 30 the room the
    # matrix keeps free at either end of the axis
    free = column - (6.2 * max(len(str(l)) for l in labels) + 30)

    return min(15, max(5, free / len(parasites)))


@st.cache_data(show_spinner=False)
def generate_tissue_dots(per_tissue, groups, group_order, niches, palette, column):
    '''
    A dot wherever a parasite is predicted to interact with a host protein expressed in a
    tissue it infects, sized by how many such interactions there are. The parasites are in
    the order of the two matrices above it -- taxonomic group, then name -- so a column is
    the same parasite everywhere on the page and a clade stays together.

    A matrix rather than the nested rectangles this used to be: the parasites of one host
    share 18 tissues between them but infect a median of two each, so what there is to see
    is which parasites meet the host in the same place -- and nesting each parasite inside
    its own rectangle is the one arrangement that never puts two of them side by side.

    The tissues are ordered by how many parasites reach them, as the host proteins of the
    dot matrix are, so the tissues every parasite has in common are the top rows.

    The parasites are named above the columns, over the same two strips the figures above
    carry. The clade is on every dot as well, but a dot is all this matrix has: a parasite
    infects a median of two tissues, so a column is two or three dots with white between
    them, and a clade read off them alone is a colour hunted for down a sparse column. The
    strip states the boundary the dots leave to be inferred.

    :param per_tissue: interactions per parasite and tissue, as
                       count_interactions_per_tissue counts them
    :param dict groups: {parasite label: taxonomic group}
    :param dict group_order: {taxonomic group: its place in the order}
    :param dict niches: {parasite label: niche}
    :param dict palette: {taxonomic group: colour}
    :param float column: pixels the figure is drawn across, which sizes the dots
    '''
    dots = per_tissue.copy()
    dots['group'] = dots['taxid1_label'].map(lambda p: groups.get(p, UNKNOWN_GROUP))
    dots['parasite'] = dots['taxid1_label'].map(lambda p: f'{p[0]}. {p.split(" ")[1]}')
    dots['niche'] = dots['taxid1_label'].map(niches).fillna(web_utils.UNKNOWN_NICHE)
    order = sorted(dots['taxid1_label'].unique(),
                   key=lambda p: (group_order.get(groups.get(p), len(group_order)), p))
    parasites = [f'{p[0]}. {p.split(" ")[1]}' for p in order]
    reach = dots.groupby('Tissue').agg(parasites=('taxid1_label', 'nunique'),
                                       total=('interactions', 'sum'))
    tissues = list(reach.sort_values(['parasites', 'total'], ascending=False, kind='stable').index)

    # the axes count cells and the names are ticks written against them, which is what lets
    # the strips be drawn: a band on a categorical axis is a category of its own, and would
    # be read as another parasite or another tissue
    dots = dots.assign(x=dots['parasite'].map({p: i for i, p in enumerate(parasites)}),
                       y=dots['Tissue'].map({t: i for i, t in enumerate(tissues)}))

    size = dot_size(parasites, tissues, column)
    sizeref = max(1, dots['interactions'].max()) / size ** 2
    figure = go.Figure()
    # one trace per group rather than one for every dot, so the groups are the legend and
    # clicking one takes that group off the plot
    for group in [g for g in list(palette) + [UNKNOWN_GROUP] if g in set(dots['group'])]:
        rows = dots[dots['group'] == group]
        figure.add_trace(go.Scatter(
            x=rows['x'], y=rows['y'], mode='markers', name=group,
            marker=dict(color=palette.get(group, UNKNOWN_COLOR), size=rows['interactions'],
                        sizemode='area', sizeref=sizeref, sizemin=4, line=dict(width=0)),
            customdata=rows[['Tissue', 'parasite', 'interactions']].to_numpy(),
            hovertemplate='%{customdata[0]}<br>%{customdata[1]}<br>predicted interactions: '
                          '%{customdata[2]}<extra></extra>'))

    groups_shown, niches_shown = add_parasite_strips(figure, dots, parasites, palette)
    height = max(420, 19 * len(tissues) + 240)
    room = add_dot_legend(figure, groups_shown, [], niches_shown, palette,
                          column - (6.2 * max(len(str(t)) for t in tissues) + 30), height)

    figure.update_layout(height=height, plot_bgcolor='white',
                         # the room the legends take is left under the plot rather than
                         # taken out of it: they are placed against the foot of the figure
                         margin=dict(l=0, r=0, t=10, b=room),
                         xaxis_title=None, yaxis_title='tissue the parasite infects')
    # every parasite is named, however narrow its column: plotly thins the labels that no
    # longer fit, and a matrix with every other column named cannot be read at all, so they
    # are drawn at the size a column has room for instead. The range runs back past the two
    # strips rather than on past the last column, which is what puts them between the names
    # and the plot
    figure.update_xaxes(range=[-0.5, len(parasites) - 0.5], side='top', tickmode='array',
                        tickvals=list(range(len(parasites))), ticktext=parasites,
                        tickangle=-60, automargin=True, ticks='',
                        tickfont=dict(size=max(SMALLEST_LABEL, min(11, round(size / 1.1)))),
                        showgrid=True, gridcolor='#f0f0f0', zeroline=False)
    # reversed, so that the tissue the most parasites reach is the top row
    figure.update_yaxes(range=[len(tissues) - 0.5, -3], tickmode='array',
                        tickvals=list(range(len(tissues))), ticktext=tissues, ticks='',
                        automargin=True, showgrid=True, gridcolor='#f0f0f0', zeroline=False)

    return figure


@st.cache_data(show_spinner=False)
def generate_cell_type_bars(per_cell_type, tissue, groups, palette):
    '''
    The cell types of one tissue, each a bar of the interactions predicted there, split by
    taxonomic group. Cell types are the level that does not fit the matrix -- 48 of them
    across the tissues, unevenly annotated -- so they are behind a choice of tissue rather
    than drawn all at once and left unreadable.

    Summing over parasites is sound here where summing over cell types is not: two
    parasites interacting in the same cell type are two different interactions.

    A bar counts the interactions with host proteins expressed above 1 nTPM in that cell
    type, which count_interactions_per_tissue defines; a protein may be counted in several,
    so the bars still overlap and do not partition the tissue.
    '''
    data = per_cell_type[per_cell_type['Tissue'] == tissue].copy()
    data['group'] = data['taxid1_label'].map(lambda p: groups.get(p, UNKNOWN_GROUP))
    data['parasite'] = data['taxid1_label'].map(lambda p: f'{p[0]}. {p.split(" ")[1]}')
    totals = data.groupby('Cell type')['interactions'].sum().sort_values(ascending=False,
                                                                        kind='stable')
    cell_types = list(totals.index)

    figure = px.bar(data, x='interactions', y='Cell type', color='group', orientation='h',
                    color_discrete_map=palette, category_orders={
                        'Cell type': cell_types,
                        'group': [g for g in palette if g in set(data['group'])]},
                    custom_data=['parasite'])
    # two parasites of one group are two segments of the same colour, and would read as a
    # single bar without a line to part them
    figure.update_traces(marker_line=dict(color='white', width=1),
                         hovertemplate='%{y}<br>%{customdata[0]}<br>predicted interactions: '
                                       '%{x}<extra></extra>')
    figure.update_layout(height=max(320, 26 * len(cell_types) + 160), plot_bgcolor='white',
                         margin=dict(l=0, r=0, t=10, b=10), legend_title_text='',
                         legend=dict(orientation='h', yanchor='bottom', y=1.01, x=0),
                         xaxis_title='predicted interactions', yaxis_title=None,
                         bargap=0.25)
    # the bars are stacked, so the axis has to reach the total of a cell type and not the
    # tallest of its parts
    web_utils.count_ticks(figure, totals.max(), showgrid=True, gridcolor='#f0f0f0')

    return figure

@st.cache_data(show_spinner=False)
def get_tissue_expressed_predictions(data_dir, config, host_taxids, score=MIN_SCORE):
    '''
    Predictions restricted to the host proteins that are expressed in a tissue the
    parasite is known to infect (config['parasites'][taxid]['tissues']), which is the
    same restriction the tissue matrix applies. What is left are the
    interactions that could take place where the parasite actually is, rather than every
    interaction predicted from orthology. Parasites left without any interactor simply
    do not appear in what is built from this.

    `score` drops the interactions predicted below that confidence, the same cut the
    network page offers. It is applied here rather than per figure so the shared-interactor
    heatmap and the shared-protein dots keep counting the same interactions.

    One row is one predicted interaction, parasite protein (`source`) included: the heatmap
    only ever looks at which host proteins are reached, but the dot matrix sizes its dots by
    how many parasite proteins reach each one.
    '''
    predictions = web_utils.get_host_predictions(data_dir, host_taxids)[['taxid1', 'taxid1_label', 'source',
                                                               'target', 'target_name', 'weight']]
    predictions = predictions[predictions['weight'] >= score].drop('weight', axis=1)
    tissues = utils.read_parquet_file(input_file=f'{data_dir}/tissues_cell_types.parquet')
    expressed = tissues.rename({'Gene': 'target'}, axis=1)[['target', 'Tissue']].drop_duplicates()

    mapped_tissues = config['tissues']
    infected_tissues = pd.DataFrame([(str(taxid), mapped_tissues[t].lower())
                                     for taxid, parasite in config['parasites'].items()
                                     for t in parasite['tissues']],
                                    columns=['taxid1', 'Tissue'])

    aux = predictions.drop_duplicates().astype({'taxid1': str})
    aux = pd.merge(aux, expressed, on='target')
    aux = pd.merge(aux, infected_tissues, on=['taxid1', 'Tissue'])

    return aux[['taxid1_label', 'source', 'target', 'target_name']].drop_duplicates()


# what the two axes of the matrix and the columns of the dot plot can be ordered on.
# Ordering is what makes a block visible: parasites next to each other are compared by
# eye, so the order decides which question the figure answers -- whether a clade shares
# its interactors, or whether the parasites inside a host cell do
ORDER_BY_GROUP = 'taxonomic group'
ORDER_BY_NICHE = 'intracellular / extracellular'


def parasite_order(labels, groups, group_order, niches, order_by):
    '''
    The parasites in the order the figures put them on their axes: the chosen annotation
    first, so its values are contiguous, then the other one, then the name. Ordering on
    the niche keeps the clades blocked inside each niche rather than scattering them, so
    the figure gains the niche blocks without losing the ones it already had.

    :param labels: the parasite labels to order
    :param dict groups: {parasite label: taxonomic group}
    :param dict group_order: {taxonomic group: its rank in the config}
    :param dict niches: {parasite label: niche}, as web_utils.get_niches builds them
    :param str order_by: ORDER_BY_GROUP or ORDER_BY_NICHE
    :return: the labels, sorted
    '''
    def rank(parasite):
        clade = group_order.get(groups.get(parasite), len(group_order))
        niche = web_utils.NICHE_ORDER.index(niches[parasite]) \
            if niches.get(parasite) in web_utils.NICHE_ORDER else len(web_utils.NICHE_ORDER)

        return (niche, clade, parasite) if order_by == ORDER_BY_NICHE else (clade, niche,
                                                                           parasite)

    return sorted(labels, key=rank)


@st.cache_data(show_spinner=False)
def get_shared_interactor_similarity(df_pred, groups, group_order, niches, order_by):
    '''
    How alike the host interactors of each pair of parasites are, as the Jaccard similarity
    of the two sets: the host proteins both parasites reach, over the host proteins either
    of them reaches. The diagonal is left empty -- a parasite is identical to itself, which
    says nothing about any pair and would be the darkest cell of every row.

    The ratio and not the count of shared proteins. A count follows how many interactors
    the two parasites have between them, so the pairs it puts at the top are the pairs of
    the best-predicted parasites: two nematodes of eight hundred interactors share more
    proteins by having more, and a pair that share almost everything they have is a pale
    cell beside them. Divided by the union, a cell is the overlap of the pair and nothing
    else, and the blocks against the diagonal are blocks of parasites reaching the same
    host proteins rather than blocks of parasites with many.

    The count is returned beside the ratio: it is what the hover of a cell gives the ratio
    in terms of, and what the dialog behind a click then lists.

    The parasites are in the order of the dot matrix, whichever annotation that order is
    taken from, so that a row is the same parasite in both and a block against the
    diagonal is a block in both.

    :return: (similarity, shared counts, clades, niches), the two frames on the same
             parasites in the same order, or None where there are fewer than three
             parasites to compare
    '''
    targets = {g: set(df['target']) for g, df in df_pred.groupby('taxid1_label')}
    targets = {g: t for g, t in targets.items() if t}
    labels = parasite_order(targets, groups, group_order, niches, order_by)
    if len(labels) < 3:
        return None

    shared = np.array([[len(targets[a] & targets[b]) for b in labels] for a in labels],
                      dtype=float)
    # the union of two non-empty sets is non-empty, and the sets are filtered to those
    union = np.array([[len(targets[a] | targets[b]) for b in labels] for a in labels],
                     dtype=float)
    similarity = shared / union
    np.fill_diagonal(similarity, np.nan)
    np.fill_diagonal(shared, np.nan)

    return (pd.DataFrame(similarity, index=labels, columns=labels),
            pd.DataFrame(shared, index=labels, columns=labels),
            [groups.get(g, UNKNOWN_GROUP) for g in labels],
            [niches.get(g, web_utils.UNKNOWN_NICHE) for g in labels])


@st.cache_data(show_spinner=False)
def flat_colour_scale(values, colors):
    '''
    The colour scale of a strip of flat bands -- one band per value, nothing interpolated
    between two of them -- and the codes that index it, which is what a strip is drawn from.

    The colour is repeated at both ends of a band so that nothing is interpolated between
    two values; the steps have to be built in order, since sorting them puts the two
    entries of a boundary in colour order rather than in band order, which hands a band its
    neighbour's colour.

    :param values: the value of each cell of the strip, in the order it is drawn
    :param dict colors: {value: colour}, in the order the bands are numbered
    :return: (heatmap keywords for the scale, the code of each value)
    '''
    steps = [step for i, value in enumerate(colors)
             for step in ([i / len(colors), colors[value]],
                          [(i + 1) / len(colors), colors[value]])]

    return dict(colorscale=steps, zmin=-0.5, zmax=len(colors) - 0.5, showscale=False,
                hovertemplate='%{text}<extra></extra>'), \
        [list(colors).index(v) for v in values]


def generate_shared_interactor_heatmap(similarity, shared, clades, niches, palette, column):
    '''
    The shared-interactor similarity matrix, with a strip of the taxonomic group and a strip
    of the niche of each parasite down the side and along the top, so that the two axes are
    visibly the same list of parasites in the same order.

    The columns are named above the matrix, outside their own strips, and the two legends
    sit below it. A name and the two colours that annotate it are then read together
    rather than from opposite sides of the matrix, and the legends, which belong to no
    part of the figure in particular, take the room nothing else wants.

    Both strips are drawn whichever of them the parasites are ordered on. The one that was
    ordered on comes out in blocks and is the question the figure is being read for; the
    other stays there to be read against it, which is where a block that does not follow
    the order is seen -- a clade whose parasites share their interactors across two niches,
    or the other way about.

    The cells are square because the figure is sized for it: both axes span the same number
    of cells and the margins leave the same number of pixels between them, so the square is
    arithmetic done here rather than a constraint plotly applies as it draws. Constraining
    the axes to their domain is the more direct way of asking, but the domain plotly settles
    on is written back into the figure, and Streamlit redraws the figure it last drew rather
    than the one this function returns: opened full screen and closed again, the matrix came
    back shrunk to the shape the full screen had wanted, and shrank again on every opening.

    The strips are drawn on the axes of the matrix rather than in subplots of their own,
    which is what keeps them against it: a subplot has a domain of its own, and would leave
    the matrix floating away from its own labels.

    The diagonal is drawn grey by a trace of its own rather than being left blank:
    blank renders as the white the colour scale starts at, so it could not be told from a
    pair sharing nothing. It carries no value, only the name of the parasite.

    Every cell off the diagonal carries an invisible marker as well, which is what a click
    on it lands on and what its hover is read from, and the index of that trace is returned
    beside the figure so the page can tell such a click from any other.

    :param similarity: the Jaccard similarity of every pair, as
                       get_shared_interactor_similarity builds it, which is what the cells
                       are coloured by
    :param shared: the host proteins each pair has in common, which the hover gives the
                   similarity of a cell in terms of
    :param column: pixels the column holding the figure is on the screen the page is being
                   read on, which web_utils.column_width measures. It is the height that is
                   sized from it -- the width belongs to the column -- so that the square
                   the cells are held to is the whole of the figure rather than a part of it
    :return: the figure, and the index of the trace holding the clickable cells
    '''
    shown = [g for g in palette if g in set(clades)]
    niches_shown = [n for n in web_utils.NICHE_ORDER + [web_utils.UNKNOWN_NICHE]
                    if n in set(niches)]
    x_names = [f'{g[0]}. {g.split(" ")[1]}' for g in similarity.index]
    y_names = list(similarity.index)
    cells = list(range(len(y_names)))
    strip, codes = flat_colour_scale(clades, {g: palette[g] for g in shown})
    niche_strip, niche_codes = flat_colour_scale(
        niches, {n: web_utils.NICHE_COLORS[n] for n in niches_shown})

    # the same span on both axes -- the matrix, the two strips beside it and the cell
    # between them and it -- so that a square plot area is one of square cells
    span = len(cells) + 2.6
    # the ticks beside the colour bar are the whole of what varies in the room it needs,
    # and a similarity is written to two decimals however large the largest of them is
    right = COLORBAR_ROOM + COLORBAR_DIGIT * len(f'{np.nanmax(similarity.to_numpy()):{TICK}}')

    # A line of text is about 1.1 times its point size tall, which is the room a cell has to
    # give the label against it. The cells being square, one size would do for both axes --
    # the names above the matrix are held to a floor of their own, being the only thing that
    # names a column and the first to be lost.
    def fits(room, smallest=6):
        sizes = [size for size in (10, 9, 8, 7, 6) if size >= smallest]
        return next(size for size in sizes if room >= 1.1 * size or size == sizes[-1])

    longest = max(len(name) for name in y_names)

    def square(size):
        '''Side of the square and the left margin, with the names written `size` points.'''
        # what the cells are held square to: what is left of the column once the labels and
        # the colour bar have taken theirs, and never larger than the rows want to be
        margin = 0.65 * size * longest + 12

        return max(240, min(CELL * span, column - margin - right)), margin

    # the room the labels need comes off the column before the square can be measured, and
    # their size follows from the square, so it is measured once at the largest the names
    # can be written and again at the size that left them
    font = fits(square(10)[0] / span)
    side, left = square(font)
    # the names above the matrix are all that names a column, so they are not written
    # smaller than they can be read: forty of them run into each other on a narrow screen
    # rather than being drawn at a size nobody can make out on any screen
    x_font = fits(side / span, smallest=SMALLEST_LABEL)

    # the room the names of the columns need above the matrix. Rotated sixty degrees they
    # stand about 0.87 of their length tall, and reach half their length out to the right of
    # the last column as well, which the margin the colour bar is given has the room for
    top = 0.87 * (0.65 * x_font * max(len(name) for name in x_names) + 12) + 25
    # room below the matrix for the two legends. Plotly wraps the entries to the width of
    # the plot, so the rows they take are counted here rather than assumed: a row that has
    # not been paid for is one plotly makes by taking it off the plot, and a plot area that
    # is no longer the square the cells are drawn in
    def legend_rows(names):
        entry = LEGEND_ENTRY + LEGEND_CHAR * max(len(name) for name in names)

        return -(-len(names) // max(1, int(side // entry)))

    clade_rows = legend_rows(shown)
    niche_rows = legend_rows(niches_shown)
    bottom = 34 + LEGEND_ROW * (clade_rows + niche_rows)

    figure = go.Figure()
    # the group and the niche of each parasite, beside its row and above its column. x0/y0
    # put a strip a cell clear of the matrix, dx/dy give it a cell of its own to fill. The
    # group is the inner strip of the two, being the one the figure has always carried
    for offset, (values, code_list, scale) in enumerate([(clades, codes, strip),
                                                         (niches, niche_codes, niche_strip)]):
        figure.add_trace(go.Heatmap(z=[[code] for code in code_list],
                                    x0=-1.4 - offset, dx=1, y=cells,
                                    text=[[v] for v in values], ygap=1, **scale))
        figure.add_trace(go.Heatmap(z=[code_list], x=cells,
                                    y0=-1.4 - offset, dy=1,
                                    text=[list(values)], xgap=1, **scale))

    # the similarities. The hover of a cell belongs to the clickable layer added below
    # rather than to this trace: two traces answering the same pointer answer it
    # differently, and the one that can be clicked is the one that should say what
    # clicking it opens.
    #
    # The scale runs to the largest similarity on the matrix rather than to the 1.0 a
    # Jaccard can reach, as it ran to the largest count before it. The pairs of a host are
    # nothing like evenly spread -- a human's run to 0.8 between sibling species, with the
    # median pair at 0.08 -- so a host whose parasites are less alike would be drawn in the
    # pale end of a fixed scale and could not be read at all. The colour bar states the
    # values either way
    figure.add_trace(go.Heatmap(z=similarity.to_numpy(), x=cells, y=cells, hoverinfo='skip',
                                colorscale=['#ffffff', '#deebf7', '#9ecae1', '#6baed6',
                                            '#3182bd', '#08519c'],
                                zmin=0, zmax=np.nanmax(similarity.to_numpy()),
                                hoverongaps=False,
                                colorbar=dict(title=dict(text='Jaccard similarity',
                                                         side='right'),
                                              thickness=12, len=0.6, y=1, yanchor='top',
                                              x=1 + COLORBAR_GAP / side,
                                              tickformat=TICK,
                                              tickfont=dict(size=10))))

    # the diagonal, a cell of flat grey per parasite, drawn over the empty cells the count
    # matrix leaves. Its hover names the parasite alone: the cell stands for no pair, so
    # there is no shared-interactor count to give. It is laid over the whole grid, one cell
    # of it drawn and the rest empty, so both traces turn hovering on gaps off -- left on,
    # the empty cells of whichever trace is on top answer for the cells beneath them
    diagonal = np.full(similarity.shape, np.nan)
    np.fill_diagonal(diagonal, 0)
    figure.add_trace(go.Heatmap(z=diagonal, x=cells, y=cells,
                                text=[[name] * len(y_names) for name in y_names],
                                colorscale=[[0, DIAGONAL_COLOUR], [1, DIAGONAL_COLOUR]],
                                zmin=0, zmax=1, showscale=False, hoverongaps=False,
                                hovertemplate='%{text}<extra></extra>'))

    # what a click lands on, and what the hover of a cell is read from. A heatmap cell
    # cannot be clicked at all: Streamlit picks a click up from the selection plotly makes
    # of it, and a heatmap is not a trace plotly can select anything in, so the click
    # reaches the page as nothing. A scatter is, so every cell of the matrix carries an
    # invisible square marker the size of the cell, which puts the whole of the cell in
    # reach of the pointer. The diagonal is left without one: it stands for no pair.
    #
    # The axes count cells rather than name parasites, so a marker is found again by the
    # cell it sits on -- its x and y are the indices of the two parasites in `counts`.
    click_targets = [(x, y) for y in cells for x in cells if x != y]
    figure.add_trace(go.Scatter(
        x=[x for x, _ in click_targets], y=[y for _, y in click_targets], mode='markers',
        marker=dict(symbol='square', size=side / span, color='rgba(0,0,0,0)',
                    line=dict(width=0)),
        # plotly dims what was not selected, which on a click would leave the one cell
        # that was clicked lit and wash the rest of the matrix out behind the dialog
        selected=dict(marker=dict(opacity=1)), unselected=dict(marker=dict(opacity=1)),
        # the count beside the ratio: a similarity of 0.12 says how alike the pair are
        # and nothing about how much there is of it, and the two of them are what the
        # dialog behind the click then lists
        text=[f'{y_names[y]} and {y_names[x]}<br>'
              f'Jaccard similarity {similarity.iat[y, x]:{TICK}}<br>'
              f'Shared interactors {shared.iat[y, x]:.0f}' for x, y in click_targets],
        hovertemplate='%{text}<extra></extra>', showlegend=False))
    click_layer = len(figure.data) - 1

    # the strips are heatmaps and cannot carry a legend of their own, so their values are
    # named by empty traces whose only purpose is their legend entry. Two legends and not
    # one: the strips are two different facts about the parasite, and run together they
    # read as one key of nine colours
    for group in shown:
        figure.add_trace(go.Scatter(x=[None], y=[None], mode='markers', name=group,
                                    marker=dict(size=10, symbol='square', color=palette[group]),
                                    hoverinfo='skip', showlegend=True))
    for niche in niches_shown:
        figure.add_trace(go.Scatter(x=[None], y=[None], mode='markers', name=niche,
                                    marker=dict(size=10, symbol='square',
                                                color=web_utils.NICHE_COLORS[niche]),
                                    legend='legend2', hoverinfo='skip', showlegend=True))

    # every parasite is named on both axes, however small the cells are drawn: plotly thins
    # tick labels that no longer fit, and a heatmap with every third row labelled cannot be
    # read at all
    ticks = dict(tickmode='array', tickvals=cells, ticks='')
    # named above the matrix, where the strips of a column are. The label starts at the
    # tick it belongs to and leans up and to the right of it, which is the slant the names
    # are drawn at under the other figures of the page: the same name is read the same way
    # wherever on the page it is met
    figure.update_xaxes(range=[-3, len(cells) - 0.4], side='top',
                        ticktext=x_names, tickangle=-60, tickfont=dict(size=x_font), **ticks)
    # reversed, so that the first parasite is the top row and the diagonal runs the way it
    # is read; the strips above the columns are drawn before the first row of the matrix,
    # so the range runs back past them rather than on past the last row
    figure.update_yaxes(range=[len(cells) - 0.4, -3],
                        ticktext=y_names, tickfont=dict(size=font), **ticks)

    # the plot area is `side` pixels each way, the margins holding the labels, the legend
    # and the colour bar out of it, and that is what makes the cells square. Plotly widens a
    # margin of its own accord where what sits in it does not fit -- a legend wrapped onto
    # more rows than there is room for below the matrix -- and the cells are then drawn a
    # little wider than they are tall, which is the whole of what a bad measurement costs
    # the niche legend sits below the clade one, a row of it clear of the plot, so the
    # rows counted into `bottom` are the rows the two of them actually take
    entries = dict(orientation='h', yanchor='top', xanchor='left', x=0, itemclick=False,
                   itemdoubleclick=False, font=dict(size=11), title_font=dict(size=11))
    figure.update_layout(width=left + side + right, height=side + top + bottom,
                         plot_bgcolor='white',
                         margin=dict(l=left, r=right, t=top, b=bottom),
                         legend=dict(y=-0.01, **entries),
                         legend2=dict(y=-0.01 - LEGEND_ROW * clade_rows / side, **entries))

    return figure, click_layer


# longest descriptive protein name written beside a gene symbol on the dot matrix. The
# label is an axis tick and every character of it is taken off the width of the plot, so
# the description is cut where it stops earning the space -- the whole of it is in the
# hover, which is where a name that long is read anyway. Forty characters take more than
# half of the column on a laptop, which leaves the forty parasites of a human a couple of
# pixels of matrix each.
DESCRIPTION_WIDTH = 26


def label_proteins(df_pred, annotations, truncate=None):
    '''
    Names each host protein by its gene symbol and its descriptive protein name, since the
    symbol on its own identifies the protein only for someone who already knows it.

    The descriptions are keyed by STRING id and the matrix is keyed by symbol, so the
    first description found for a symbol is the one used. Proteins UniProt has nothing
    but "Uncharacterized protein" for keep their symbol alone: that description names
    nothing and would be repeated down the axis.

    :param df_pred: tissue-expressed predictions of the host group
    :param dict annotations: STRING id --> descriptive protein name
    :param int truncate: characters of the description a label has room for, or None for
                         the whole of it
    :return: {gene symbol: protein label}
    '''
    described = {}
    for protein, name in df_pred[['target', 'target_name']].drop_duplicates().values:
        description = (annotations or {}).get(protein, '')
        if description and not description.lower().startswith('uncharacterized'):
            described.setdefault(name, description.rstrip('.'))

    labels = {}
    for name in df_pred['target_name'].unique():
        description = described.get(name)
        if truncate and description and len(description) > truncate:
            description = description[:truncate - 1].rstrip() + '…'
        labels[name] = f'{name} · {description}' if description else str(name)

    return labels


def summarise_localisations(df_pred, localisations):
    '''
    The DeepLoc call of each host protein of the matrix, keyed by the gene symbol the
    rows are drawn under rather than by the STRING id DeepLoc wrote it against.

    A host protein reaches this page only because DeepLoc put it somewhere a parasite of
    its host can reach, but those places are not the same interaction: a cell membrane
    protein is met on the surface of the host cell, an extracellular one in the fluid
    around it, and a cytosolic or nuclear one only by a parasite with an intracellular
    stage, which is inside the cell to meet it. Which of them a protein was kept for is
    web_utils.classify_localisation, which the home page reads both sides with.

    Read on all four classes, and not on the ones the niche of a parasite allows, as the
    home page reads a host protein. A row of the figures here is a protein across every
    parasite that reaches it, of either niche, so there is no one niche to read it on: what
    the band says is where DeepLoc puts the protein, and the strip above the columns is
    where the parasite is that meets it there.

    A gene can have more than one STRING protein identifier, so the id DeepLoc is most
    sure of is the one the row is described by.

    :param df_pred: tissue-expressed predictions of the host group
    :param localisations: DeepLoc table written by pipeline/build_deeploc_localisations.py
    :return: dataframe indexed by gene symbol, or None if there are no localisations
    '''
    if localisations is None or localisations.empty:
        return None

    proteins = df_pred[['target', 'target_name']].drop_duplicates()
    called = proteins.merge(localisations, left_on='target', right_on='protein', how='inner')
    if called.empty:
        return None

    scores = [c for c in web_utils.DEEPLOC_SCORES.values() if c in called]
    called['best'] = called[scores].max(axis=1)
    called = called.sort_values('best', ascending=False, kind='stable')
    called = called.drop_duplicates('target_name').set_index('target_name')

    called['surface'] = web_utils.classify_localisation(called)

    return called[['surface'] + scores + ['localizations']]


@st.dialog('Host interactors shared by a pair of parasites', width='large')
def show_shared_interactors_dialog(first, second, df_pred, annotations, localisations):
    '''
    The host proteins behind one cell of the shared-interactor matrix: the cell gives a
    count, and this is what the count is made of.

    Opened over the matrix rather than placed under it, as the network page opens the
    AlphaFold models of an interaction, so that the figures below do not move on every
    click.

    The rows are host proteins and not gene symbols, since that is what the cell counts:
    a gene with more than one STRING identifier is more than one interactor of the
    matrix, and collapsing the two would put a number in the dialog that the cell the
    dialog was opened from disagrees with. The identifier is written beside the name,
    which is where the two rows of such a gene are told apart.

    Each parasite has a column of its own, holding the number of its proteins predicted to
    reach that host protein -- being reached by eighty-six proteins of a parasite and by
    one are not the same prediction, and the pair is the whole subject here.

    :param first: parasite of the row of the cell that was clicked
    :param second: parasite of its column
    :param df_pred: tissue-expressed predictions, as the matrix counted them, so the rows
                    answer to the number in the cell
    :param dict annotations: STRING id --> descriptive protein name
    :param localisations: DeepLoc table, which the surface class is read from
    '''
    edges = df_pred[df_pred['taxid1_label'].isin([first, second])]
    targets = {p: set(df['target']) for p, df in edges.groupby('taxid1_label')}
    shared = targets.get(first, set()) & targets.get(second, set())

    st.markdown(f'**{first}** and **{second}**', unsafe_allow_html=True)
    st.caption(f'The {len(shared)} host proteins both parasites are predicted to interact '
               f'with, at the confidence the page is set to and among the proteins '
               f'expressed in a tissue each parasite infects. A column counts the proteins '
               f'of that parasite predicted to reach the host protein.')
    if not shared:
        return

    edges = edges[edges['target'].isin(shared)]
    # one column per parasite, indexed by the host protein: both parasites reach every
    # protein here, the rows being the ones they share
    degree = edges.groupby(['taxid1_label', 'target'])['source'].nunique().unstack('taxid1_label')
    names = edges.drop_duplicates('target').set_index('target')['target_name']
    labels = label_proteins(edges, annotations)

    table = pd.DataFrame({'Host protein': [labels.get(names[t], names[t]) for t in degree.index],
                          'Identifier': list(degree.index)})
    columns = []
    for parasite in (first, second):
        # named as the axes of the matrix name a parasite, so a column is read back to the
        # row or the column of the cell the dialog was opened from
        column_name = f'{parasite[0]}. {parasite.split(" ")[1]} proteins'
        table[column_name] = degree[parasite].fillna(0).astype(int).values
        columns.append(column_name)

    surface = summarise_localisations(edges, localisations)
    if surface is not None:
        table['DeepLoc'] = [surface['surface'].get(names[t], NO_LOCALISATION)
                            for t in degree.index]

    # the proteins reached by the most parasite proteins first: the rows the pair has most
    # of are the rows the pair is being read for
    table = table.assign(_reach=table[columns].sum(axis=1)).sort_values(
        ['_reach', 'Host protein'], ascending=[False, True], kind='stable').drop(columns='_reach')

    st.dataframe(table, width='stretch', hide_index=True)
    st.download_button('Download table', table.to_csv(index=False).encode('utf-8'),
                       file_name=f'shared_interactors_{first}_{second}.csv'.replace(' ', '_'),
                       mime='text/csv')


@st.cache_data(show_spinner=False)
def get_top_shared_proteins(df_pred, groups, group_order, niches, order_by,
                            annotations=None, localisations=None, top=40):
    '''
    The host proteins that the most parasites are predicted to interact with. The heatmap
    counts how much two parasites have in common but does not say what they have in common,
    which is what this is for. Proteins only one parasite interacts
    with are left out: they are not shared by anything.

    Each dot also carries `degree`, the number of proteins of that parasite predicted to
    interact with that host protein -- the degree of the host protein in the network of
    that one parasite. A dot is otherwise only a yes, and one parasite protein reaching a
    host protein is a thinner prediction than eighty-six of them.

    :param dict annotations: STRING id --> descriptive protein name, which the rows are
                             labelled with beside the gene symbol
    :param localisations: DeepLoc table, which the dots carry the surface class and the
                          two surface probabilities of the host protein from
    :return: the dots, the row labels, the column labels, and how many host proteins more
             than one parasite reaches -- the rows are the `top` of those, and the caption
             says so where that leaves some of them out
    '''
    edges = df_pred[['taxid1_label', 'source', 'target', 'target_name']].drop_duplicates()
    # One row per dot, keyed by the display gene name rather than a protein identifier, so
    # aliases do not produce dots on top of each other.
    pairs = edges[['taxid1_label', 'target_name']].drop_duplicates()
    counts = pairs.groupby('target_name')['taxid1_label'].nunique()
    counts = counts[counts > 1].sort_values(ascending=False, kind='stable')
    if counts.empty:
        return None

    proteins = list(counts.head(top).index)
    degree = edges.groupby(['taxid1_label', 'target_name'])['source'].nunique()
    dots = pairs[pairs['target_name'].isin(proteins)].copy()
    dots['group'] = dots['taxid1_label'].map(lambda p: groups.get(p, UNKNOWN_GROUP))
    dots['parasites'] = dots['target_name'].map(counts)
    dots['degree'] = pd.MultiIndex.from_frame(dots[['taxid1_label', 'target_name']]).map(degree)
    dots['parasite'] = dots['taxid1_label'].map(lambda p: f'{p[0]}. {p.split(" ")[1]}')
    labels = label_proteins(df_pred, annotations, DESCRIPTION_WIDTH)
    dots['protein'] = dots['target_name'].map(labels)
    full_labels = label_proteins(df_pred, annotations)
    dots['protein_full'] = dots['target_name'].map(full_labels)
    surface = summarise_localisations(df_pred, localisations)
    if surface is not None:
        dots = dots.join(surface, on='target_name')
        dots['surface'] = dots['surface'].fillna(NO_LOCALISATION)
        dots['localizations'] = dots['localizations'].fillna('')
    order = parasite_order(dots['taxid1_label'].unique(), groups, group_order, niches,
                           order_by)
    dots['niche'] = dots['taxid1_label'].map(niches).fillna(web_utils.UNKNOWN_NICHE)

    return (dots, [labels[p] for p in proteins],
            [f'{p[0]}. {p.split(" ")[1]}' for p in order], len(counts))


def localisation_labels(dots, proteins):
    '''
    What the band beside a row of the dot plot says when it is hovered: the class of the
    host protein, and for a protein called for more than one, which classes those are.

    `Several` is one colour over two different statements -- a protein on the cell membrane
    and in the fluid around it, and one in the cytosol and the nucleus, which are opposite
    sides of that membrane -- and a legend of one entry cannot tell them apart. The band
    names them instead, where a reader who wonders what a purple cell is made of asks.

    :param dots: the frame behind the figure, carrying `surface` and the probability of
                 each class beside each protein
    :param list proteins: the host proteins in the order of the rows
    :return: the label of each row, in that order
    '''
    called = dots.drop_duplicates('protein').set_index('protein')
    scores = {c: web_utils.DEEPLOC_SCORES[c] for c in web_utils.HOST_CLASSES
              if web_utils.DEEPLOC_SCORES[c] in called}
    labels = []
    for protein in proteins:
        row = called.loc[protein]
        crossed = [c for c, column in scores.items()
                   if row[column] > web_utils.DEEPLOC_CUTOFFS[c]]
        # a protein of one class is named by it, and one the filter had nothing to say
        # about -- no localisation at all -- by what the legend calls it
        labels.append(' + '.join(crossed) if len(crossed) > 1 else row['surface'])

    return labels


def add_parasite_strips(figure, dots, parasites, palette):
    '''
    The taxonomic group and the niche of each parasite, as two bands above its column.
    Drawn on both dot matrices of the page and on the shared-interactor heatmap, so that
    the three figures the parasites run across are read the same way round and a column is
    found by the same two colours wherever it is met.

    Drawn as heatmap cells inside the axes rather than as shapes hung off them. A shape has
    to be placed against the plot area, and the plot area is whatever the row labels leave
    of the column once automargin has taken the room they need; a cell is placed on the
    axes themselves, at negative coordinates the range is opened up to hold, and follows
    the plot wherever the labels leave it.

    The group is the inner strip of the two, as on the heatmap: it is the one every figure
    has always carried, and the two are read in the same order on all of them.

    :param figure: the figure, modified in place
    :param dots: the frame behind it, carrying `group` and `niche` beside each `parasite`
    :param list parasites: the parasites in the order of the x axis
    :param dict palette: {taxonomic group: colour}
    :return: (groups, niches) drawn, each in the order it is keyed in
    '''
    clade_of = dict(zip(dots['parasite'], dots['group']))
    niche_of = dict(zip(dots['parasite'], dots['niche']))
    clades = [clade_of.get(p, UNKNOWN_GROUP) for p in parasites]
    niches = [niche_of.get(p, web_utils.UNKNOWN_NICHE) for p in parasites]
    groups_shown = [g for g in list(palette) + [UNKNOWN_GROUP] if g in set(clades)]
    niches_shown = [n for n in web_utils.NICHE_ORDER + [web_utils.UNKNOWN_NICHE]
                    if n in set(niches)]

    # the group a cell clear of the first row and the niche outside it
    for offset, (values, colors) in enumerate(
            [(clades, {g: palette.get(g, UNKNOWN_COLOR) for g in groups_shown}),
             (niches, {n: web_utils.NICHE_COLORS[n] for n in niches_shown})]):
        scale, codes = flat_colour_scale(values, colors)
        figure.add_trace(go.Heatmap(z=[codes], x=list(range(len(parasites))),
                                    y0=-1.4 - offset, dy=1, text=[values], xgap=1, **scale))

    return groups_shown, niches_shown


def add_localisation_strip(figure, dots, proteins):
    '''
    Where DeepLoc puts each host protein, as a band beside its row. It is a property of the
    protein and so constant along the row, which is what makes it a band rather than a
    channel of the dots: everything about a parasite is above the plot and everything about
    a protein beside it, each axis annotated by the thing it lists.

    :param figure: the dot plot, modified in place
    :param dots: the frame behind it, carrying `surface` beside each `protein`
    :param list proteins: the host proteins in the order of the y axis
    :return: the localization classes drawn, in the order they are keyed in
    '''
    surface_of = dict(zip(dots['protein'], dots['surface']))
    localisations = [surface_of.get(p, NO_LOCALISATION) for p in proteins]
    shown = [c for c in LOCALISATION_ORDER if c in set(localisations)]
    scale, codes = flat_colour_scale(localisations,
                                     {c: LOCALISATION_COLORS[c] for c in shown})
    # the colour of a cell is the class, its hover the classes behind that class
    figure.add_trace(go.Heatmap(z=[[code] for code in codes], x0=-1.4, dx=1,
                                y=list(range(len(proteins))),
                                text=[[label] for label in
                                      localisation_labels(dots, proteins)],
                                ygap=1, **scale))

    return shown


def add_dot_legend(figure, groups, localisations, niches, palette, width, height):
    '''
    Names the two strips the dot plot carries no legend for, and lays the three keys of the
    figure out under it as three legends of their own.

    The taxonomic groups are the colour of the dots themselves and are named by the traces
    that draw them -- clicking one still takes that group off the plot -- while a strip is a
    heatmap, which carries no legend entry of its own, so its values are named by empty
    traces whose only purpose is the entry they leave behind. A square, being what a band of
    colour is keyed by.

    Three legends and not one. The three channels have values that read alike -- a parasite
    outside the host cell and a host protein outside it are both "Extracellular" -- and run
    together as one key of fourteen colours there is nothing to say which of them a colour
    belongs to. Named and set apart, each key is read against the strip it annotates.

    The room each legend takes is counted rather than assumed, since the next legend is
    placed under the last: plotly wraps the entries to the width of the plot, and a row
    that has not been paid for is a row drawn over by the legend below it. The title is
    paid for as a row of its own, being where plotly puts it on a horizontal legend, and
    the entries wrap inside the whole width rather than what a title beside them would
    leave.

    :param figure: the dot plot, modified in place
    :param list groups: the taxonomic groups drawn, which the dot traces already name
    :param list localisations: the localization classes drawn, in the order of the legend
    :param list niches: the niches drawn, in the order of the legend
    :param dict palette: {taxonomic group: colour}
    :param float width: pixels of the plot the legends wrap inside
    :param float height: pixels the figure is drawn down, which the offsets are a share of
    :return: the pixels the three legends take, which is the room to leave under the plot
    '''
    # the dots of a group are sized by their degree, so the entry plotly would write for
    # such a trace carries the size of whichever dot came first -- a Cestoda of degree one
    # is a speck beside a Nematoda of eighty-six. The traces are taken out of the legend and
    # named by a marker of one size instead, keeping their legendgroup so that clicking the
    # entry still takes the group off the plot
    for trace in figure.data:
        trace.update(legendgroup=trace.name, showlegend=False)
    for name, color in [(g, palette.get(g, UNKNOWN_COLOR)) for g in groups]:
        figure.add_scatter(x=[None], y=[None], mode='markers', name=name, legendgroup=name,
                           marker=dict(symbol='circle', size=10, color=color),
                           hoverinfo='skip', showlegend=True)
    for name, color, legend in ([(c, LOCALISATION_COLORS[c], 'legend2') for c in localisations]
                                + [(n, web_utils.NICHE_COLORS[n], 'legend3') for n in niches]):
        figure.add_scatter(x=[None], y=[None], mode='markers', name=name, legend=legend,
                           marker=dict(symbol='square', size=10, color=color),
                           hoverinfo='skip', showlegend=True)

    def pixels(names):
        '''The room a legend takes: its title, and the rows its entries wrap onto.'''
        entry = LEGEND_ENTRY + DOT_LEGEND_CHAR * max(len(name) for name in names)

        return LEGEND_TITLE_ROW + LEGEND_ROW * -(-len(names)
                                                 // max(1, int(width // entry)))

    keys = [(legend, title, names, pixels(names))
            for legend, title, names in [('legend', 'taxonomic group', groups),
                                         ('legend2', 'DeepLoc', localisations),
                                         ('legend3', web_utils.NICHE_TITLE, niches)]
            if names]
    room = sum(taken for *_, taken in keys) + LEGEND_PAD

    # measured from the foot of the figure and not from the foot of the plot: a legend
    # placed against the plot moves with the margin the rotated names take off the top,
    # which is measured by plotly as it draws and is not known here, and three legends that
    # move by a margin nobody counted are three legends drawn over each other
    offset = room / height
    for legend, title, names, taken in keys:
        figure.update_layout({legend: dict(
            orientation='h', yanchor='top', y=offset, yref='container', x=0,
            traceorder='normal', title_text=title, title_font=dict(size=11),
            font=dict(size=11))})
        offset -= taken / height

    return room


@st.cache_data(show_spinner=False)
def generate_shared_protein_dots(dots, proteins, parasites, palette, column):
    '''
    A dot wherever a parasite is predicted to interact with one of the proteins, the
    parasites in the order the other figures use so the taxonomic groups stay together, and
    the proteins ordered by how many parasites reach them. The dot is sized by how many
    proteins of that parasite reach that host protein.

    The area of the dot is what carries the degree (as plotly express drew it), since that
    is the channel size is read on, and sizemin keeps the single-protein dots -- the
    largest group of them -- from collapsing to a speck next to a degree of eighty-six.

    The dot is left with the one shape and the one channel. Where DeepLoc puts the host
    protein runs beside the row instead, in the colours the home page splits its bars in:
    it is a property of the protein, so as a marker shape it was drawn once for every
    parasite reaching it, at sizes where a circle and a diamond are the same dot, and it
    took the area of that marker away from the degree it is supposed to say.

    The parasites are named above the columns, with the group and the niche between the
    names and the plot, which is how the matrix beside this figure is laid out: the two
    figures are the same parasites in the same order, and are read the same way round.

    The axes count cells rather than name a parasite or a protein, the names being ticks
    written against them. Categorical axes would place a dot as readily, but the strips
    could not be drawn on them: a band on a categorical axis is a category of its own, and
    would be read as another parasite or another protein.

    :param dots: one row per predicted (parasite, host protein) pair, as
                 get_top_shared_proteins builds them
    :param list proteins: the host proteins, most-shared first, as the rows are labelled
    :param list parasites: the parasites in the order of the columns
    :param dict palette: {taxonomic group: colour}
    :param float column: pixels the figure is drawn across, which sizes the dots
    '''
    localised = 'surface' in dots.columns
    cells = {'x': {p: i for i, p in enumerate(parasites)},
             'y': {p: i for i, p in enumerate(proteins)}}
    dots = dots.assign(x=dots['parasite'].map(cells['x']),
                       y=dots['protein'].map(cells['y'])).dropna(subset=['x', 'y'])

    # the hover is written out rather than left to a column name: the protein and the
    # parasite are the two axes already, and `parasites` and `degree` are two different
    # counts of two different things, which as bare numbers under their column names they
    # do not say. The axes count cells, so the parasite is named from the row behind the
    # dot rather than from the x it sits on
    hover_columns = ['protein_full', 'parasites', 'degree', 'parasite']
    hover_lines = ['%{customdata[0]}', 'parasites reaching it: %{customdata[1]}',
                   'proteins of %{customdata[3]} reaching it: %{customdata[2]}']
    if localised:
        # the class is the band beside the row, so the hover carries the probabilities
        # behind it rather than naming it a second time: the four classes the host filter
        # reads, two to a line, the surface of the cell first and the inside of it second
        scores = [web_utils.DEEPLOC_SCORES[c] for c in web_utils.HOST_CLASSES
                  if web_utils.DEEPLOC_SCORES[c] in dots]
        parts = [f'{DEEPLOC_LABELS[name]} %{{customdata[{len(hover_columns) + i}]:.2f}}'
                 for i, name in enumerate(scores)]
        hover_columns += scores
        hover_lines += [', '.join(parts[i:i + 2]) for i in range(0, len(parts), 2)]

    size = dot_size(parasites, proteins, column)
    # the area of the largest dot stands for the largest degree, which is the size plotly
    # express solved for when it drew this figure
    sizeref = max(1, dots['degree'].max()) / size ** 2
    figure = go.Figure()
    # one trace per group rather than one for every dot, so the groups are the legend and
    # clicking one takes that group off the plot
    for group in [g for g in list(palette) + [UNKNOWN_GROUP] if g in set(dots['group'])]:
        rows = dots[dots['group'] == group]
        figure.add_trace(go.Scatter(
            x=rows['x'], y=rows['y'], mode='markers', name=group,
            marker=dict(color=palette.get(group, UNKNOWN_COLOR), size=rows['degree'],
                        sizemode='area', sizeref=sizeref, sizemin=4, line=dict(width=0)),
            customdata=rows[hover_columns].to_numpy(),
            hovertemplate='<br>'.join(hover_lines) + '<extra></extra>'))

    groups_shown, niches = add_parasite_strips(figure, dots, parasites, palette)
    localisations = add_localisation_strip(figure, dots, proteins) if localised else []
    # what the legends wrap inside is what the row labels leave of the column, which is what
    # the dots were sized on as well
    height = max(420, 19 * len(proteins) + 240)
    room = add_dot_legend(figure, groups_shown, localisations, niches, palette,
                          column - (6.2 * max(len(str(p)) for p in proteins) + 30), height)

    figure.update_layout(height=height, plot_bgcolor='white',
                         # the room the legends take is left under the plot rather than
                         # taken out of it: they are placed against the foot of the figure
                         margin=dict(l=0, r=0, t=10, b=room),
                         xaxis_title=None, yaxis_title='host protein')
    # every parasite is named, however narrow its column: plotly thins the labels that no
    # longer fit, and a matrix with every other column named cannot be read at all, so they
    # are drawn at the size a column has room for instead. The names lean up and to the
    # right of the tick they belong to, the slant they are drawn at everywhere on the page.
    # The range runs back past the two strips above the columns rather than on past the
    # last of them, which is what puts them between the names and the plot
    figure.update_xaxes(range=[-1.9 if localised else -0.5, len(parasites) - 0.5],
                        side='top', tickmode='array', tickvals=list(range(len(parasites))),
                        ticktext=parasites, tickangle=-60, automargin=True, ticks='',
                        tickfont=dict(size=max(SMALLEST_LABEL, min(11, round(size / 1.1)))),
                        showgrid=True, gridcolor='#f0f0f0', zeroline=False)
    # reversed, so that the most-shared protein is the top row. The labels are as long as a
    # protein name, so plotly is left to take the room they need off the plot rather than
    # drawing them over it or cutting them at the margin
    figure.update_yaxes(range=[len(proteins) - 0.5, -3], tickmode='array',
                        tickvals=list(range(len(proteins))), ticktext=proteins, ticks='',
                        automargin=True, showgrid=True, gridcolor='#f0f0f0', zeroline=False)

    return figure


# The figures below sit two to a row and are stretched to their column, which is all most
# of them need. The shared-interactor matrix has to know how wide that column came out --
# its cells are square, so its height is its width -- so the browser is asked as the page
# opens and answers on the run after: until then the figures are drawn for a laptop
# (web_utils.DEFAULT_PAGE_WIDTH) and the matrix is a square with room to spare above it.
column = web_utils.column_width(2)

st.caption('The parasites predicted against one host, compared with each other: which host '
           'interactors they share, which host proteins several of them reach, and the '
           'tissues and cell types in which their interactions can take place.')
st.markdown("---")

col1, col2, col3 = st.columns(3)

with col1:
    st.write('')

with col2:
    # The host selection is shared across pages.
    selected_host, selected_taxids = web_utils.host_selector(
        config, web_utils.load_predictions(data_dir),
        'Select a host to compare the parasites that infect it')
    if selected_host == web_utils.NO_HOST:
        st.text('Choose 1 host to explore the predicted host-parasite interactions')

with col3:
    st.write('')


if selected_host != web_utils.NO_HOST:
    # The heatmap and dot matrix read the same filtered interactions: the heatmap gives
    # every parasite pair a shared-host-interactor count, while the dot matrix names the
    # host proteins those pairs have in common.
    parasite_groups = {p['label']: p.get('group', UNKNOWN_GROUP)
                       for p in config['parasites'].values()}
    group_order = {g: i for i, g in enumerate(config.get('parasite_groups', {}))}
    niches = web_utils.get_niches(config)

    # one slider for the three figures below it, which are the same interactions read three
    # ways -- filtering them apart would put different numbers in figures whose captions
    # say they match. The tissue dots at the foot of the page are a different count and
    # keep every prediction.
    # It sits in the middle of three columns, as the host selector above it does: left to
    # itself a slider takes the whole width of the page, which is a metre of track for a
    # range of half a point
    slider_column, order_column = st.columns([2, 1])
    with slider_column:
        score = st.slider('Confidence score', MIN_SCORE, MAX_SCORE, DEFAULT_SCORE,
                          help='Interactions predicted below this confidence are left out of '
                               'the three figures below. The tissue plot at the foot of the '
                               'page counts every prediction.')
    with order_column:
        # ordering the parasites is what decides which blocks the two figures can show:
        # the values of whichever annotation is chosen come out contiguous, so a set of
        # parasites sharing their interactors is a square against the diagonal rather than
        # a scatter of cells to be found by reading the labels
        order_by = st.radio('Order the parasites by', [ORDER_BY_GROUP, ORDER_BY_NICHE],
                            horizontal=True,
                            help='Which annotation the two figures below put next to each '
                                 'other on their axes. The strips beside the axes show both '
                                 'either way; this is which of them comes out in blocks.')
    counted = get_tissue_expressed_predictions(data_dir, config, selected_taxids, score)
    shared_similarity = get_shared_interactor_similarity(counted, parasite_groups,
                                                        group_order, niches, order_by)
    top_shared = get_top_shared_proteins(counted, parasite_groups, group_order, niches,
                                         order_by,
                                         web_utils.load_protein_annotations(data_dir),
                                         web_utils.load_deeploc_localisations(data_dir))

    matrix, shared = st.columns(2)

    with matrix:
        st.subheader("Overlap of the host interactors of each pair of parasites")
        st.caption('How alike the host interactors of each pair of parasites are, as the '
                   'Jaccard similarity of the two sets: the host proteins both reach, over '
                   'the host proteins either of them reaches. A count of shared proteins on '
                   'its own follows how many interactors the pair have between them, and '
                   'ranks the best-predicted parasites above the most alike; this does not. '
                   'Two strips '
                   'run along each axis: the taxonomic group of the parasite, and whether it '
                   'lives inside a host cell or outside one. Whichever the parasites are '
                   'ordered by comes out in blocks. The diagonal, where a parasite meets '
                   'itself, is greyed out. Hover a cell for the count behind its ratio, or '
                   'click it to see the shared host proteins.')
        if shared_similarity is not None:
            # the figure is given the width of the column and keeps its cells square within
            # it, so it follows whatever screen the page is read on
            # the figure carries the width of the column rather than being stretched to it,
            # which is what keeps a square matrix square: stretched, the plot area is squared
            # by plotly on the fly, and it does not undo that when the figure is given its
            # column back after being opened full screen
            figure, click_layer = generate_shared_interactor_heatmap(
                *shared_similarity, config.get('parasite_groups', {}), column)
            # The chart is remounted after every dialog, its key carrying a counter. Two
            # things would otherwise keep the same cell from being opened twice running:
            # Streamlit drops a selection identical to the one it is already holding, and
            # plotly reads a second click on a selected point as a deselection. A key that
            # has changed is a chart holding no selection at all, so the next click on any
            # cell is a new one.
            nonce = st.session_state.get(CELL_NONCE_KEY, 0)
            clicked = st.plotly_chart(
                figure, width='content', on_select='rerun', selection_mode='points',
                key=f'shared_cells_{selected_host}_{order_by}_{score}_{nonce}')
            # only the cells can be clicked, but the figure carries traces that could grow
            # points of their own, so the layer the click came from is checked
            cell = next((point for point
                         in (clicked or {}).get('selection', {}).get('points', [])
                         if point.get('curve_number') == click_layer), None)
            if cell is not None:
                # the axes of the matrix count cells, so the two parasites are the row and
                # the column the marker sits on
                parasites = list(shared_similarity[0].index)
                st.session_state[CELL_NONCE_KEY] = nonce + 1
                show_shared_interactors_dialog(
                    parasites[int(cell['y'])], parasites[int(cell['x'])], counted,
                    web_utils.load_protein_annotations(data_dir),
                    web_utils.load_deeploc_localisations(data_dir))
        else:
            st.text(f'Fewer than three parasites of {selected_host} share any host protein')

    with shared:
        if top_shared is not None:
            # the figure draws the rows it is given; the count of what they were taken from
            # belongs to the caption, which is the only place saying what is on screen and
            # what is not
            *figure_arguments, shareable = top_shared
            shown = len(figure_arguments[1])
            st.subheader("Host interactors common to several parasites")
            # a truncated figure says what it is a top of: with hundreds of proteins tied
            # a few parasites apart, a reader who is not told 40 of 810 reads the last row
            # as the last protein several parasites reach
            selection = (f'The {shown} host proteins reached by the most parasites, of the '
                         f'{shareable} reached by more than one.' if shareable > shown else
                         f'The {shown} host proteins reached by more than one parasite, '
                         'most first.')
            st.caption(selection + ' A dot wherever a parasite is predicted to interact with '
                       "one, sized by the number of that parasite's proteins reaching it. "
                       'The band beside each row gives the DeepLoc 2 localization '
                       'of the host protein; '
                       'hover a dot for the probabilities behind it. Above the columns run the '
                       'taxonomic group of the parasite and whether it lives inside a host cell '
                       'or outside one.')
            st.plotly_chart(
                generate_shared_protein_dots(*figure_arguments,
                                             config.get('parasite_groups', {}), column),
                width='stretch')

    per_tissue, per_cell_type = count_interactions_per_tissue(data_dir, config,
                                                              selected_taxids, score)
    ranked = per_tissue.groupby('Tissue')['interactions'].sum().sort_values(ascending=False,
                                                                           kind='stable')
    annotated = set(per_cell_type['Tissue'])
    choices = [t for t in ranked.index if t in annotated]
    tissues, cell_types = st.columns(2)
    with tissues:
        st.subheader("Tissues in which the predicted interactions can take place")
        st.caption('Predicted interactions per parasite and tissue, restricted to the tissues '
                   'each parasite is known to infect and sized by the number of interactions with '
                   'proteins expressed there. Tissues are ordered by the number of parasites '
                   'infecting them. Above the columns run the taxonomic group of the parasite and '
                   'whether it lives inside a host cell or outside one, the same two strips the '
                   'figures above carry. An interaction is counted once per tissue, irrespective '
                   'of the number of cell types the host protein is expressed in.')
        if per_tissue.empty:
            st.info('No predicted interaction is left at this confidence in a tissue the '
                    'parasites are known to infect.')
        else:
            st.plotly_chart(generate_tissue_dots(per_tissue, parasite_groups, group_order,
                                                 niches, config.get('parasite_groups', {}),
                                                 column),
                            width='stretch')

    with cell_types:
        if choices:
            st.subheader("Cell types of a tissue")
            st.caption('Predicted interactions per cell type of the selected tissue, stacked by '
                       'taxonomic group. A host protein counts towards a cell type where its '
                       'expression exceeds 1 nTPM. A protein expressed above that threshold '
                       'in several cell types counts in each, so the bars are not a partition '
                       'of the tissue.')
            tissue = st.selectbox('Tissue', choices, index=0,
                                  help='Tissues with cell type annotation, most interactions first')
            st.plotly_chart(generate_cell_type_bars(per_cell_type, tissue, parasite_groups,
                                                    config.get('parasite_groups', {})),
                            width='stretch')

st.markdown("---")



# Footer
with st.container():
    web_utils.footer()
