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
# the height of one row of a legend above that heatmap, in pixels, which is what the room
# left above the matrix is counted in
LEGEND_ROW = 22
# the room one entry of a legend takes across the top of that heatmap, in pixels:
# its marker and the padding around it, and about seven and a half pixels a character of the
# longest name in it -- plotly gives every entry the room the longest of them needs
LEGEND_ENTRY = 40
LEGEND_CHAR = 7.5
# the smallest the parasite names under a matrix are written, in points
SMALLEST_LABEL = 9


# marker each surface class is drawn with in the shared-interactors matrix, and the order
# they are offered in the legend. Colour there is already the parasite's taxonomic group
# and size is the degree, so the localisation of the host protein goes on the shape --
# constant down a row, since it is a property of the protein and not of the interaction.
# NO_LOCALISATION is a protein DeepLoc was never run on, or one whose data directory
# predates pipeline/build_deeploc_localisations.py
NO_LOCALISATION = 'Not available'
# All of them filled: the dots are drawn without an outline, which is what an open
# symbol is made of
SURFACE_SYMBOLS = {web_utils.CELL_MEMBRANE: 'circle', web_utils.EXTRACELLULAR: 'diamond',
                   web_utils.BOTH_SURFACE: 'hexagon', web_utils.NOT_SURFACE: 'square',
                   NO_LOCALISATION: 'cross'}
# grey the surface classes are drawn in the legend in. The markers on the plot carry the
# colour of the parasite's group, so a colour here would say the class had one
SURFACE_LEGEND_COLOR = '#525252'
# names the DeepLoc columns are read under in the hover of the shared-interactors matrix
DEEPLOC_LABELS = {'surface': 'DeepLoc', 'cell_membrane': 'P(cell membrane)',
                  'extracellular': 'P(extracellular)', 'localizations': 'localizations'}


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
def generate_tissue_dots(per_tissue, groups, group_order, palette, column):
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
    '''
    dots = per_tissue.copy()
    dots['group'] = dots['taxid1_label'].map(lambda p: groups.get(p, UNKNOWN_GROUP))
    dots['parasite'] = dots['taxid1_label'].map(lambda p: f'{p[0]}. {p.split(" ")[1]}')
    order = sorted(dots['taxid1_label'].unique(),
                   key=lambda p: (group_order.get(groups.get(p), len(group_order)), p))
    parasites = [f'{p[0]}. {p.split(" ")[1]}' for p in order]
    reach = dots.groupby('Tissue').agg(parasites=('taxid1_label', 'nunique'),
                                       total=('interactions', 'sum'))
    tissues = list(reach.sort_values(['parasites', 'total'], ascending=False, kind='stable').index)

    size = dot_size(parasites, tissues, column)
    figure = px.scatter(dots, x='parasite', y='Tissue', color='group',
                        # plotly express flips category_orders on a y axis, so most-reached
                        # first puts the tissue the most parasites infect in the top row
                        color_discrete_map=palette, category_orders={
                            'parasite': parasites, 'Tissue': tissues,
                            'group': [g for g in palette if g in set(dots['group'])]},
                        size='interactions', size_max=size,
                        custom_data=['interactions'])
    figure.update_traces(marker=dict(sizemin=4, line=dict(width=0)),
                         hovertemplate='%{y}<br>%{x}<br>predicted interactions: '
                                       '%{customdata[0]}<extra></extra>')
    figure.update_layout(height=max(420, 19 * len(tissues) + 240), plot_bgcolor='white',
                         margin=dict(l=0, r=0, t=10, b=10), legend_title_text='',
                         legend=dict(orientation='h', yanchor='bottom', y=1.01, x=0),
                         xaxis_title=None, yaxis_title='tissue the parasite infects')
    # every parasite is named, however narrow its column: plotly thins the labels that no
    # longer fit, and a matrix with every other column named cannot be read at all, so they
    # are drawn at the size a column has room for instead
    figure.update_xaxes(tickangle=-60, tickmode='linear', dtick=1, automargin=True,
                        tickfont=dict(size=max(SMALLEST_LABEL, min(11, round(size / 1.1)))),
                        showgrid=True, gridcolor='#f0f0f0')
    figure.update_yaxes(automargin=True, showgrid=True, gridcolor='#f0f0f0')

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
def get_shared_interactor_counts(df_pred, groups, group_order, niches, order_by):
    '''
    Number of host proteins each pair of parasites is predicted to interact with. The
    diagonal is left empty -- a parasite shares nothing with itself that is worth counting
    against a pair -- and every off-diagonal cell is the host proteins shared by that pair.

    The parasites are in the order of the dot matrix, whichever annotation that order is
    taken from, so that a row is the same parasite in both and a block against the
    diagonal is a block in both.
    '''
    targets = {g: set(df['target']) for g, df in df_pred.groupby('taxid1_label')}
    targets = {g: t for g, t in targets.items() if t}
    labels = parasite_order(targets, groups, group_order, niches, order_by)
    if len(labels) < 3:
        return None

    counts = np.array([[len(targets[a] & targets[b]) for b in labels] for a in labels],
                      dtype=float)
    np.fill_diagonal(counts, np.nan)

    return (pd.DataFrame(counts, index=labels, columns=labels),
            [groups.get(g, UNKNOWN_GROUP) for g in labels],
            [niches.get(g, web_utils.UNKNOWN_NICHE) for g in labels])


@st.cache_data(show_spinner=False)
def generate_shared_interactor_heatmap(counts, clades, niches, palette, column):
    '''
    The shared-interactor count matrix, with a strip of the taxonomic group and a strip of
    the niche of each parasite down the side and along the bottom, so that the two axes are
    visibly the same list of parasites in the same order.

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
    pair sharing nothing. It carries no count, only the name of the parasite.

    :param column: pixels the column holding the figure is on the screen the page is being
                   read on, which web_utils.column_width measures. It is the height that is
                   sized from it -- the width belongs to the column -- so that the square
                   the cells are held to is the whole of the figure rather than a part of it
    '''
    def flat_scale(values, colors):
        '''
        One flat band per value, and the codes that index it. The colour is repeated at
        both ends of a band so that nothing is interpolated between two values; the steps
        have to be built in order, since sorting them puts the two entries of a boundary in
        colour order rather than in band order, which hands a band its neighbour's colour.
        '''
        steps = [step for i, value in enumerate(colors)
                 for step in ([i / len(colors), colors[value]],
                              [(i + 1) / len(colors), colors[value]])]

        return dict(colorscale=steps, zmin=-0.5, zmax=len(colors) - 0.5, showscale=False,
                    hovertemplate='%{text}<extra></extra>'), \
            [list(colors).index(v) for v in values]

    shown = [g for g in palette if g in set(clades)]
    niches_shown = [n for n in web_utils.NICHE_ORDER + [web_utils.UNKNOWN_NICHE]
                    if n in set(niches)]
    x_names = [f'{g[0]}. {g.split(" ")[1]}' for g in counts.index]
    y_names = list(counts.index)
    cells = list(range(len(y_names)))
    strip, codes = flat_scale(clades, {g: palette[g] for g in shown})
    niche_strip, niche_codes = flat_scale(
        niches, {n: web_utils.NICHE_COLORS[n] for n in niches_shown})

    # the same span on both axes -- the matrix, the two strips beside it and the cell
    # between them and it -- so that a square plot area is one of square cells
    span = len(cells) + 2.6
    # the counts beside the colour bar are the whole of what varies in the room it needs
    right = COLORBAR_ROOM + COLORBAR_DIGIT * len(f'{np.nanmax(counts.to_numpy()):.0f}')

    # A line of text is about 1.1 times its point size tall, which is the room a cell has to
    # give the label against it. The cells being square, one size would do for both axes --
    # the names below the matrix are held to a floor of their own, being the only thing that
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
    # the names below the matrix are all that names a column, so they are not written
    # smaller than they can be read: forty of them run into each other on a narrow screen
    # rather than being drawn at a size nobody can make out on any screen
    x_font = fits(side / span, smallest=SMALLEST_LABEL)

    bottom = 0.87 * (0.65 * x_font * max(len(name) for name in x_names) + 12) + 25
    # room above the matrix for the two legends. Plotly wraps the entries to the width of
    # the plot, so the rows they take are counted here rather than assumed: a row that has
    # not been paid for is one plotly makes by taking it off the plot, and a plot area that
    # is no longer the square the cells are drawn in
    def legend_rows(names):
        entry = LEGEND_ENTRY + LEGEND_CHAR * max(len(name) for name in names)

        return -(-len(names) // max(1, int(side // entry)))

    clade_rows = legend_rows(shown)
    niche_rows = legend_rows(niches_shown)
    top = 34 + LEGEND_ROW * (clade_rows + niche_rows)

    figure = go.Figure()
    # the group and the niche of each parasite, beside its row and under its column. x0/y0
    # put a strip a cell clear of the matrix, dx/dy give it a cell of its own to fill. The
    # group is the inner strip of the two, being the one the figure has always carried
    for offset, (values, code_list, scale) in enumerate([(clades, codes, strip),
                                                         (niches, niche_codes, niche_strip)]):
        figure.add_trace(go.Heatmap(z=[[code] for code in code_list],
                                    x0=-1.4 - offset, dx=1, y=cells,
                                    text=[[v] for v in values], ygap=1, **scale))
        figure.add_trace(go.Heatmap(z=[code_list], x=cells,
                                    y0=len(cells) + 0.4 + offset, dy=1,
                                    text=[list(values)], xgap=1, **scale))

    # the axes count cells rather than name parasites, the strips having to sit a cell out
    # from the matrix, so the pair a cell stands for is carried in its hover text
    figure.add_trace(go.Heatmap(z=counts.to_numpy(), x=cells, y=cells,
                                text=[[f'{row} and {other}' for other in y_names]
                                      for row in y_names],
                                colorscale=['#ffffff', '#deebf7', '#9ecae1', '#6baed6',
                                            '#3182bd', '#08519c'],
                                zmin=0, zmax=np.nanmax(counts.to_numpy()),
                                hovertemplate='%{text}<br>Shared interactors %{z:.0f}'
                                              '<extra></extra>', hoverongaps=False,
                                colorbar=dict(title=dict(text='Shared interactors',
                                                         side='right'),
                                              thickness=12, len=0.6, y=1, yanchor='top',
                                              x=1 + COLORBAR_GAP / side,
                                              tickfont=dict(size=10))))

    # the diagonal, a cell of flat grey per parasite, drawn over the empty cells the count
    # matrix leaves. Its hover names the parasite alone: the cell stands for no pair, so
    # there is no shared-interactor count to give. It is laid over the whole grid, one cell
    # of it drawn and the rest empty, so both traces turn hovering on gaps off -- left on,
    # the empty cells of whichever trace is on top answer for the cells beneath them
    diagonal = np.full(counts.shape, np.nan)
    np.fill_diagonal(diagonal, 0)
    figure.add_trace(go.Heatmap(z=diagonal, x=cells, y=cells,
                                text=[[name] * len(y_names) for name in y_names],
                                colorscale=[[0, DIAGONAL_COLOUR], [1, DIAGONAL_COLOUR]],
                                zmin=0, zmax=1, showscale=False, hoverongaps=False,
                                hovertemplate='%{text}<extra></extra>'))

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
    figure.update_xaxes(range=[-3, len(cells) - 0.4],
                        ticktext=x_names, tickangle=-60, tickfont=dict(size=x_font), **ticks)
    # reversed, so that the first parasite is the top row and the diagonal runs the way it
    # is read; the strips along the bottom are the last rows of the range, not the first
    figure.update_yaxes(range=[len(cells) + 2, -0.6],
                        ticktext=y_names, tickfont=dict(size=font), **ticks)

    # the plot area is `side` pixels each way, the margins holding the labels, the legend
    # and the colour bar out of it, and that is what makes the cells square. Plotly widens a
    # margin of its own accord where what sits in it does not fit -- a legend wrapped onto
    # more rows than there is room for above the matrix -- and the cells are then drawn a
    # little wider than they are tall, which is the whole of what a bad measurement costs
    # the niche legend sits above the clade one, a row of it clear of the plot, so the
    # rows counted into `top` are the rows the two of them actually take
    entries = dict(orientation='h', yanchor='bottom', xanchor='left', x=0, itemclick=False,
                   itemdoubleclick=False, font=dict(size=11), title_font=dict(size=11))
    figure.update_layout(width=left + side + right, height=side + top + bottom,
                         plot_bgcolor='white',
                         margin=dict(l=left, r=right, t=top, b=bottom),
                         legend=dict(y=1.01, **entries),
                         legend2=dict(y=1.01 + LEGEND_ROW * clade_rows / side, **entries))

    return figure


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

    A host protein reaches this page only because DeepLoc called it surface-exposed, but
    the ways of being surface-exposed are not the same interaction: a cell membrane protein
    is met on the surface of the host cell, an extracellular one is met in the fluid around
    it, and a protein DeepLoc assigns both classes is met either way. Which of them a
    protein was kept for is web_utils.classify_surface, which the home page reads the
    parasite proteins with.

    A gene can have more than one STRING protein identifier, so the id DeepLoc is most
    sure is surface-exposed is the one the row is described by.

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

    called['best'] = called[['cell_membrane', 'extracellular']].max(axis=1)
    called = called.sort_values('best', ascending=False, kind='stable')
    called = called.drop_duplicates('target_name').set_index('target_name')

    called['surface'] = web_utils.classify_surface(called)

    return called[['surface', 'cell_membrane', 'extracellular', 'localizations']]


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
            [f'{p[0]}. {p.split(" ")[1]}' for p in order], int(counts.max()))


def split_dot_legend(figure, groups, surfaces, palette, niches=()):
    '''
    Splits the legend of the shared-interactors matrix back into the channels it is
    drawn from.

    Plotly express draws a trace per combination of the two channels it is given and names
    it after both, so a matrix coloured by taxonomic group and shaped by surface class
    comes out with a legend of "Nematoda, Cell membrane" entries, one per combination that
    occurs. Every drawn trace is taken out of the legend and the legend is built from
    entries carrying no data instead: one per taxonomic group, in its colour and always as
    a circle so the shape of whichever combination came first says nothing, one per
    surface class, in grey since the shape means the same whatever the colour it is drawn
    in, and one per niche of the strip above the columns. The group entries keep the
    legendgroup of the traces they name, so clicking one still hides that parasite group.

    :param figure: the matrix figure, modified in place
    :param list groups: taxonomic groups that occur, in the order of the legend
    :param list surfaces: surface classes that occur, in the order of the legend
    :param dict palette: {taxonomic group: colour}
    :param niches: niches that occur, in the order of the legend, named by a square: the
                   strip above the columns is a band of solid colour, and a square is what
                   a band of colour is keyed by. Its greys are the greys of the strip and
                   not the one grey the surface classes are drawn in, the colour being the
                   whole of what a niche is told by
    '''
    for trace in figure.data:
        trace.update(legendgroup=trace.name.split(',')[0].strip(), showlegend=False)

    def add_key(name, symbol, color, legendgroup):
        figure.add_scatter(x=[None], y=[None], mode='markers', name=name,
                           marker=dict(symbol=symbol, size=9, color=color),
                           legendgroup=legendgroup, hoverinfo='skip', showlegend=True)

    for group in groups:
        add_key(group, 'circle', palette.get(group, UNKNOWN_COLOR), group)
    for surface in surfaces:
        add_key(surface, SURFACE_SYMBOLS[surface], SURFACE_LEGEND_COLOR, 'surface')
    for niche in niches:
        add_key(niche, 'square', web_utils.NICHE_COLORS[niche], 'niche')


# where the niche strip of the dot plot sits, as a fraction of the height of the plot:
# its own height, and the gap between it and the top of the plot. Shallow enough to read
# as a strip about the columns rather than as a row of the matrix, which is what a band of
# the height of a row of dots would be taken for
STRIP_HEIGHT = 0.022
STRIP_GAP = 0.006


def add_niche_strip(figure, dots, parasites):
    '''
    The niche of each parasite as a band of colour across the top of the dot plot, drawn as
    shapes rather than as a trace.

    A trace would have to sit on the axes of the plot, and both of them are categorical and
    already spoken for -- the host proteins down one and the parasites across the other --
    so a band drawn on them would be a row of the matrix and read as one. Shapes hang off
    the axes instead, on the columns in x and on the plot itself in y.

    Above the columns and not under them, which is where the matrix beside it carries its
    own. The names of the parasites are under the columns, rotated and as long as a species
    name; a strip below them is a strip adrift from what it annotates, and one above them
    is a strip drawn over them.

    :param figure: the dot plot, modified in place
    :param dots: the frame behind it, carrying a `niche` beside each `parasite`
    :param list parasites: the parasites in the order of the x axis
    :return: the niches drawn, in the order they are read in
    '''
    niche_of = dict(zip(dots['parasite'], dots['niche']))
    for position, parasite in enumerate(parasites):
        niche = niche_of.get(parasite, web_utils.UNKNOWN_NICHE)
        # the columns are categories, so a band is placed by the index of its column and
        # spans the half cell either side of it that the column itself takes
        figure.add_shape(type='rect', xref='x', yref='paper', layer='below',
                         x0=position - 0.5, x1=position + 0.5,
                         y0=1 + STRIP_GAP, y1=1 + STRIP_GAP + STRIP_HEIGHT,
                         fillcolor=web_utils.NICHE_COLORS[niche], line_width=0)

    return [n for n in web_utils.NICHE_ORDER + [web_utils.UNKNOWN_NICHE]
            if n in set(niche_of.values())]


@st.cache_data(show_spinner=False)
def generate_shared_protein_dots(dots, proteins, parasites, most, palette, column):
    '''
    A dot wherever a parasite is predicted to interact with one of the proteins, the
    parasites in the order the other figures use so the taxonomic groups stay together, and the
    proteins ordered by how many parasites reach them. The dot is sized by how many
    proteins of that parasite reach that host protein.

    The area of the dot is what carries the degree (plotly's default), since that is the
    channel size is read on, and sizemin keeps the single-protein dots -- the largest group
    of them -- from collapsing to a speck next to a degree of eighty-six.

    Where the dots carry a DeepLoc call the shape of the marker is the surface class of
    the host protein, which is constant down a row: a circle is met on the cell membrane,
    a diamond in the fluid around the cell, and a hexagon either way, DeepLoc having
    assigned it both. The hover carries the two probabilities behind that and everywhere
    else DeepLoc puts the protein.

    Under the columns runs the same niche strip the matrix beside it carries, so that the
    two figures are read as the same parasites in the same order. The colour of a dot is
    already spent on the clade and its shape on the localization of the host protein, so
    the niche has nothing left to be drawn in and goes under the columns instead.
    '''
    localised = 'surface' in dots.columns
    orders = {'parasite': parasites, 'protein': proteins,
              'group': [g for g in palette if g in set(dots['group'])]}
    # the hover is written out rather than left to plotly express, which prints the raw
    # column name of whatever it is given. The protein and the parasite are the two axes
    # already, and `parasites` and `degree` are two different counts of two different
    # things, which as bare numbers under their column names they do not say
    hover_columns = ['protein_full', 'parasites', 'degree']
    hover_lines = ['%{customdata[0]}', 'parasites reaching it: %{customdata[1]}',
                   'proteins of %{x} reaching it: %{customdata[2]}']
    if localised:
        orders['surface'] = [s for s in SURFACE_SYMBOLS if s in set(dots['surface'])]
        hover_columns += ['cell_membrane', 'extracellular']
        # the class is the shape of the dot, so the hover carries the two probabilities
        # behind it rather than naming it a second time
        hover_lines.append('P(cell membrane) %{customdata[3]:.2f}, '
                           'P(extracellular) %{customdata[4]:.2f}')

    size = dot_size(parasites, proteins, column)
    figure = px.scatter(dots, x='parasite', y='protein', color='group',
                        # plotly express flips category_orders on a y axis, so `proteins`
                        # most-shared first puts the most-shared protein in the top row
                        color_discrete_map=palette, category_orders=orders,
                        symbol='surface' if localised else None,
                        symbol_map=SURFACE_SYMBOLS if localised else {},
                        size='degree', size_max=size, labels=DEEPLOC_LABELS,
                        custom_data=hover_columns)
    figure.update_traces(marker=dict(sizemin=4, line=dict(width=0)),
                         hovertemplate='<br>'.join(hover_lines) + '<extra></extra>')
    # the strip goes on before the legend is built: split_dot_legend takes every drawn
    # trace out of the legend, and the keys naming the strip are drawn traces of their own
    niches_shown = add_niche_strip(figure, dots, parasites)
    if localised:
        split_dot_legend(figure, orders['group'], orders['surface'], palette, niches_shown)
    else:
        # nothing to split the legend apart from, so the clade entries plotly express drew
        # are left as they are and the strip is named beside them
        for niche in niches_shown:
            figure.add_scatter(x=[None], y=[None], mode='markers', name=niche,
                               marker=dict(symbol='square', size=9,
                                           color=web_utils.NICHE_COLORS[niche]),
                               legendgroup='niche', hoverinfo='skip', showlegend=True)
    figure.update_layout(height=max(420, 19 * len(proteins) + 240), plot_bgcolor='white',
                         margin=dict(l=0, r=0, t=10, b=10), legend_title_text='',
                         # clear of the strip above the columns, which the legend would
                         # otherwise be drawn across
                         legend=dict(orientation='h', yanchor='bottom',
                                     y=1 + 2 * STRIP_GAP + STRIP_HEIGHT, x=0),
                         xaxis_title=None, yaxis_title=f'host protein (up to {most} parasites)')
    # the labels are as long as a protein name, so plotly is left to take the room they
    # need off the plot rather than drawing them over it or cutting them at the margin
    figure.update_yaxes(automargin=True)
    # every parasite is named, however narrow its column: plotly thins the labels that no
    # longer fit, and a matrix with every other column named cannot be read at all, so they
    # are drawn at the size a column has room for instead. They are rotated, so they are as
    # tall as they are long and are cut at the foot of the figure unless automargin takes
    # the room they need off the plot
    figure.update_xaxes(tickangle=-60, tickmode='linear', dtick=1, automargin=True,
                        tickfont=dict(size=max(SMALLEST_LABEL, min(11, round(size / 1.1)))),
                        showgrid=True, gridcolor='#f0f0f0')
    figure.update_yaxes(showgrid=True, gridcolor='#f0f0f0')

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
    shared_counts = get_shared_interactor_counts(counted, parasite_groups, group_order,
                                                niches, order_by)
    top_shared = get_top_shared_proteins(counted, parasite_groups, group_order, niches,
                                         order_by,
                                         web_utils.load_protein_annotations(data_dir),
                                         web_utils.load_deeploc_localisations(data_dir))

    matrix, shared = st.columns(2)

    with matrix:
        st.subheader("Host interactors shared by each pair of parasites")
        st.caption('Number of host interactors shared by each pair of parasites. Two strips '
                   'run along each axis: the taxonomic group of the parasite, and whether it '
                   'lives inside a host cell or outside one. Whichever the parasites are '
                   'ordered by comes out in blocks. The diagonal, where a parasite meets '
                   'itself, is greyed out.')
        if shared_counts is not None:
            # the figure is given the width of the column and keeps its cells square within
            # it, so it follows whatever screen the page is read on
            # the figure carries the width of the column rather than being stretched to it,
            # which is what keeps a square matrix square: stretched, the plot area is squared
            # by plotly on the fly, and it does not undo that when the figure is given its
            # column back after being opened full screen
            st.plotly_chart(generate_shared_interactor_heatmap(
                *shared_counts, config.get('parasite_groups', {}), column), width='content')
        else:
            st.text(f'Fewer than three parasites of {selected_host} share any host protein')

    with shared:
        if top_shared is not None:
            st.subheader("Host interactors common to several parasites")
            st.caption('Host proteins reached by the most parasites, with a dot wherever a '
                       'parasite is predicted to interact with one, sized by the number of that '
                       "parasite's proteins reaching it. Proteins reached by a single parasite "
                       'are omitted. Dot shape gives the DeepLoc 2 localization of the host '
                       'protein, circles cell membrane, diamonds extracellular, and hexagons '
                       'both; hover for the underlying probabilities. The strip under the '
                       'columns shows whether the parasite lives inside a host cell or '
                       'outside one.')
            st.plotly_chart(
                generate_shared_protein_dots(*top_shared, config.get('parasite_groups', {}),
                                             column),
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
                   'infecting them. An interaction is counted once per tissue, irrespective of '
                   'the number of cell types the host protein is expressed in.')
        if per_tissue.empty:
            st.info('No predicted interaction is left at this confidence in a tissue the '
                    'parasites are known to infect.')
        else:
            st.plotly_chart(generate_tissue_dots(per_tissue, parasite_groups, group_order,
                                                 config.get('parasite_groups', {}), column),
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
