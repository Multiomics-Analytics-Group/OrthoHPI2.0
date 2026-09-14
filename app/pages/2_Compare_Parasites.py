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

config = utils.read_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()


# fallback for a parasite without a `group` in the config
UNKNOWN_GROUP = 'Unclassified'
UNKNOWN_COLOR = '#999999'
# same default and range as the network page
MIN_SCORE, MAX_SCORE, DEFAULT_SCORE = 0.35, 0.9, 0.35
# side of a heatmap cell in pixels, which sizes the shared-interactor figure
CELL = 22
# the diagonal of the heatmap: a grey, so it is not read as a count
DIAGONAL_COLOUR = '#9e9e9e'
# pixels the colour bar takes right of the matrix; held in pixels so plotly does not widen
# the margin and unsquare the plot
COLORBAR_GAP = 8
COLORBAR_ROOM = 71
COLORBAR_DIGIT = 6.5
# two decimals: most pairs sit under a tenth
TICK = '.2f'
# height of one legend row below the heatmap, in pixels
LEGEND_ROW = 22
# room one legend entry takes, in pixels: marker, padding and ~7.5 px a character of the
# longest name
LEGEND_ENTRY = 40
LEGEND_CHAR = 7.5
# room left under the last legend of a dot matrix, in pixels
LEGEND_PAD = 14
# ~6 px a character of the longest name in a dot-matrix legend
DOT_LEGEND_CHAR = 6.0
# room a horizontal legend's title takes above its entries
LEGEND_TITLE_ROW = 18
# smallest the parasite names on a matrix axis are written, in points
SMALLEST_LABEL = 9
# pixels per row of the shared-interactor dot plot
DOT_ROW = 19
DOT_CHROME = 240
# session key counting remounts of the shared-interactor matrix, so the same cell can be
# opened twice running
CELL_NONCE_KEY = 'shared_cell_nonce'
# the tissue selector value meaning no tissue filter
ALL_TISSUES = 'All tissues'


# localization of a host protein, drawn as a band beside its row; NO_LOCALISATION is a
# protein DeepLoc was never run on
NO_LOCALISATION = 'Not available'
LOCALISATION_ORDER = list(web_utils.HOST_CLASSES) + [web_utils.SEVERAL,
                                                     web_utils.NOT_SURFACE, NO_LOCALISATION]
# paler than the grey of a rejected protein
LOCALISATION_COLORS = {**web_utils.LOCALISATION_COLORS, NO_LOCALISATION: '#f0f0f0'}
# names the DeepLoc columns are read under in the hover
DEEPLOC_LABELS = {'surface': 'DeepLoc', 'localizations': 'localizations',
                  **{column: f'P({name.lower()})'
                     for name, column in web_utils.DEEPLOC_SCORES.items()}}


@st.cache_data(show_spinner=False)
def count_interactions_per_tissue(data_dir, config, host_taxids, score=MIN_SCORE):
    '''
    Predicted interactions per parasite and tissue, and per parasite, tissue and cell
    type, keeping only the tissues each parasite is known to infect
    (config['parasites']).
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

    # keyed by display gene name so aliases of the same host gene do not inflate the counts
    def pairs(frame, *level):
        return (frame.drop_duplicates(list(level) + ['source', 'target_name'])
                     .groupby(list(level), observed=True).size()
                     .rename('interactions').reset_index())

    return (pairs(aux, 'taxid1_label', 'Tissue'),
            pairs(web_utils.keep_expressed_cell_types(aux), 'taxid1_label', 'Tissue',
                  'Cell type'))


def plot_room(labels, column):
    '''
    Pixels of a `column`-wide figure left for the plot itself, once rows named by
    `labels` have taken the room their names need.
    '''
    # 6.2 px a character of the default axis font, 30 px free at either end of the axis
    return column - (6.2 * max(len(str(l)) for l in labels) + 30)


def dot_size(parasites, labels, column):
    '''
    Diameter of the largest dot of a matrix of `parasites` columns whose rows are named
    by `labels`, drawn in a column `column` pixels wide.
    '''
    return min(15, max(5, plot_room(labels, column) / len(parasites)))


@st.cache_data(show_spinner=False)
def generate_tissue_dots(per_tissue, groups, group_order, niches, palette, column):
    '''
    A dot wherever a parasite is predicted to interact with a host protein expressed in
    a tissue it infects, sized by how many such interactions there are.
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

    # the axes count cells and the names are ticks, so the strips can be drawn without
    # becoming a category
    dots = dots.assign(x=dots['parasite'].map({p: i for i, p in enumerate(parasites)}),
                       y=dots['Tissue'].map({t: i for i, t in enumerate(tissues)}))

    size = dot_size(parasites, tissues, column)
    sizeref = max(1, dots['interactions'].max()) / size ** 2
    figure = go.Figure()
    # one trace per group, so the groups are the legend
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
                         # the legends are placed against the foot of the figure, under the
                         # plot
                         margin=dict(l=0, r=0, t=10, b=room),
                         xaxis_title=None, yaxis_title='tissue the parasite infects')
    # every parasite is named (plotly thins labels that no longer fit); the range runs back
    # past the two strips
    figure.update_xaxes(range=[-0.5, len(parasites) - 0.5], side='top', tickmode='array',
                        tickvals=list(range(len(parasites))), ticktext=parasites,
                        tickangle=-60, automargin=True, ticks='',
                        tickfont=dict(size=max(SMALLEST_LABEL, min(11, round(size / 1.1)))),
                        showgrid=True, gridcolor='#f0f0f0', zeroline=False)
    # reversed, so the tissue the most parasites reach is the top row
    figure.update_yaxes(range=[len(tissues) - 0.5, -3], tickmode='array',
                        tickvals=list(range(len(tissues))), ticktext=tissues, ticks='',
                        automargin=True, showgrid=True, gridcolor='#f0f0f0', zeroline=False)

    return figure


@st.cache_data(show_spinner=False)
def generate_cell_type_bars(per_cell_type, tissue, groups, palette):
    '''
    The cell types of one tissue, each a bar of the interactions predicted there, split
    by taxonomic group.
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
    # two parasites of one group would read as a single bar without a line to part them
    figure.update_traces(marker_line=dict(color='white', width=1),
                         hovertemplate='%{y}<br>%{customdata[0]}<br>predicted interactions: '
                                       '%{x}<extra></extra>')
    figure.update_layout(height=max(320, 26 * len(cell_types) + 160), plot_bgcolor='white',
                         margin=dict(l=0, r=0, t=10, b=10), legend_title_text='',
                         legend=dict(orientation='h', yanchor='bottom', y=1.01, x=0),
                         xaxis_title='predicted interactions', yaxis_title=None,
                         bargap=0.25)
    # the bars are stacked, so the axis has to reach the total
    web_utils.count_ticks(figure, totals.max(), showgrid=True, gridcolor='#f0f0f0')

    return figure

@st.cache_data(show_spinner=False)
def get_tissue_expressed_predictions(data_dir, config, host_taxids, score=MIN_SCORE,
                                     tissue=None):
    '''
    Predictions restricted to the host proteins that are expressed in a tissue the
    parasite is known to infect (config['parasites'][taxid]['tissues']), which is the
    same restriction the tissue matrix applies.
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
    if tissue is not None:
        aux = aux[aux['Tissue'] == tissue]

    return aux[['taxid1_label', 'source', 'target', 'target_name']].drop_duplicates()


# what the matrix axes and dot plot columns can be ordered on; the order decides which
# blocks are visible
ORDER_BY_GROUP = 'taxonomic group'
ORDER_BY_NICHE = 'intracellular / extracellular'


def parasite_order(labels, groups, group_order, niches, order_by):
    '''
    The parasites in the order the figures put them on their axes: the chosen annotation
    first, so its values are contiguous, then the other one, then the name.
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
    How alike the host interactors of each pair of parasites are, as the Jaccard
    similarity of the two sets: the host proteins both parasites reach, over the host
    proteins either of them reaches.
    '''
    targets = {g: set(df['target']) for g, df in df_pred.groupby('taxid1_label')}
    targets = {g: t for g, t in targets.items() if t}
    labels = parasite_order(targets, groups, group_order, niches, order_by)
    if len(labels) < 3:
        return None

    shared = np.array([[len(targets[a] & targets[b]) for b in labels] for a in labels],
                      dtype=float)
    # the union of two non-empty sets is non-empty
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
    The colour scale of a strip of flat bands -- one band per value, nothing
    interpolated between two of them -- and the codes that index it, which is what a
    strip is drawn from.
    '''
    steps = [step for i, value in enumerate(colors)
             for step in ([i / len(colors), colors[value]],
                          [(i + 1) / len(colors), colors[value]])]

    return dict(colorscale=steps, zmin=-0.5, zmax=len(colors) - 0.5, showscale=False,
                hovertemplate='%{text}<extra></extra>'), \
        [list(colors).index(v) for v in values]


def generate_shared_interactor_heatmap(similarity, shared, clades, niches, palette, column):
    '''
    The shared-interactor similarity matrix, with a strip of the taxonomic group and a
    strip of the niche of each parasite down the side and along the top, so that the two
    axes are visibly the same list of parasites in the same order.
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

    # the same span on both axes, so a square plot area is one of square cells
    span = len(cells) + 2.6
    # the ticks beside the colour bar are all that varies in the room it needs
    right = COLORBAR_ROOM + COLORBAR_DIGIT * len(f'{np.nanmax(similarity.to_numpy()):{TICK}}')

    # a line of text is about 1.1 times its point size tall
    def fits(room, smallest=6):
        sizes = [size for size in (10, 9, 8, 7, 6) if size >= smallest]
        return next(size for size in sizes if room >= 1.1 * size or size == sizes[-1])

    longest = max(len(name) for name in y_names)

    def square(size):
        '''Side of the square and the left margin, with the names written `size` points.'''
        # what the cells are held square to
        margin = 0.65 * size * longest + 12

        return max(240, min(CELL * span, column - margin - right)), margin

    # measured once at the largest size and again at the size that left them
    font = fits(square(10)[0] / span)
    side, left = square(font)
    # column names are held to a floor of their own, being the only thing that names a
    # column
    x_font = fits(side / span, smallest=SMALLEST_LABEL)

    # rotated sixty degrees the names stand about 0.87 of their length tall
    top = 0.87 * (0.65 * x_font * max(len(name) for name in x_names) + 12) + 25
    # plotly wraps legend entries to the plot width, so the rows are counted rather than
    # assumed
    def legend_rows(names):
        entry = LEGEND_ENTRY + LEGEND_CHAR * max(len(name) for name in names)

        return -(-len(names) // max(1, int(side // entry)))

    clade_rows = legend_rows(shown)
    niche_rows = legend_rows(niches_shown)
    bottom = 34 + LEGEND_ROW * (clade_rows + niche_rows)

    figure = go.Figure()
    # group and niche strips beside the rows and above the columns, a cell clear of the
    # matrix
    for offset, (values, code_list, scale) in enumerate([(clades, codes, strip),
                                                         (niches, niche_codes, niche_strip)]):
        figure.add_trace(go.Heatmap(z=[[code] for code in code_list],
                                    x0=-1.4 - offset, dx=1, y=cells,
                                    text=[[v] for v in values], ygap=1, **scale))
        figure.add_trace(go.Heatmap(z=[code_list], x=cells,
                                    y0=-1.4 - offset, dy=1,
                                    text=[list(values)], xgap=1, **scale))

    # hover belongs to the clickable layer below. The scale runs to the largest similarity
    # rather than 1.0
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

    # the diagonal, drawn over the empty cells; both traces turn hovering on gaps off so
    # neither answers for the cells beneath
    diagonal = np.full(similarity.shape, np.nan)
    np.fill_diagonal(diagonal, 0)
    figure.add_trace(go.Heatmap(z=diagonal, x=cells, y=cells,
                                text=[[name] * len(y_names) for name in y_names],
                                colorscale=[[0, DIAGONAL_COLOUR], [1, DIAGONAL_COLOUR]],
                                zmin=0, zmax=1, showscale=False, hoverongaps=False,
                                hovertemplate='%{text}<extra></extra>'))

    # a heatmap cannot be clicked in Streamlit, so every cell carries an invisible scatter
    # marker; x and y are the indices of the two parasites
    click_targets = [(x, y) for y in cells for x in cells if x != y]
    figure.add_trace(go.Scatter(
        x=[x for x, _ in click_targets], y=[y for _, y in click_targets], mode='markers',
        marker=dict(symbol='square', size=side / span, color='rgba(0,0,0,0)',
                    line=dict(width=0)),
        # plotly would otherwise dim everything but the clicked cell
        selected=dict(marker=dict(opacity=1)), unselected=dict(marker=dict(opacity=1)),
        # the count beside the ratio, which is what the dialog then lists
        text=[f'{y_names[y]} and {y_names[x]}<br>'
              f'Jaccard similarity {similarity.iat[y, x]:{TICK}}<br>'
              f'Shared interactors {shared.iat[y, x]:.0f}' for x, y in click_targets],
        hovertemplate='%{text}<extra></extra>', showlegend=False))
    click_layer = len(figure.data) - 1

    # the strips are heatmaps and cannot carry a legend, so their values are named by empty
    # traces
    for group in shown:
        figure.add_trace(go.Scatter(x=[None], y=[None], mode='markers', name=group,
                                    marker=dict(size=10, symbol='square', color=palette[group]),
                                    hoverinfo='skip', showlegend=True))
    for niche in niches_shown:
        figure.add_trace(go.Scatter(x=[None], y=[None], mode='markers', name=niche,
                                    marker=dict(size=10, symbol='square',
                                                color=web_utils.NICHE_COLORS[niche]),
                                    legend='legend2', hoverinfo='skip', showlegend=True))

    # every parasite is named, however small the cells
    ticks = dict(tickmode='array', tickvals=cells, ticks='')
    # named above the matrix, slanted the way the other figures of the page slant them
    figure.update_xaxes(range=[-3, len(cells) - 0.4], side='top',
                        ticktext=x_names, tickangle=-60, tickfont=dict(size=x_font), **ticks)
    # reversed, so the first parasite is the top row; the range runs back past the strips
    figure.update_yaxes(range=[len(cells) - 0.4, -3],
                        ticktext=y_names, tickfont=dict(size=font), **ticks)

    # the plot area is `side` pixels each way, the margins holding everything else, which
    # keeps the cells square
    entries = dict(orientation='h', yanchor='top', xanchor='left', x=0, itemclick=False,
                   itemdoubleclick=False, font=dict(size=11), title_font=dict(size=11))
    figure.update_layout(width=left + side + right, height=side + top + bottom,
                         plot_bgcolor='white',
                         margin=dict(l=left, r=right, t=top, b=bottom),
                         legend=dict(y=-0.01, **entries),
                         legend2=dict(y=-0.01 - LEGEND_ROW * clade_rows / side, **entries))

    return figure, click_layer


# longest description written beside a gene symbol on the dot matrix; the whole of it is in
# the hover
DESCRIPTION_WIDTH = 26


def label_proteins(df_pred, annotations, truncate=None):
    '''
    Names each host protein by its gene symbol and its descriptive protein name, since
    the symbol on its own identifies the protein only for someone who already knows it.
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
    # one column per parasite, indexed by the shared host protein
    degree = edges.groupby(['taxid1_label', 'target'])['source'].nunique().unstack('taxid1_label')
    names = edges.drop_duplicates('target').set_index('target')['target_name']
    labels = label_proteins(edges, annotations)

    table = pd.DataFrame({'Host protein': [labels.get(names[t], names[t]) for t in degree.index],
                          'Identifier': list(degree.index)})
    columns = []
    for parasite in (first, second):
        # named as the matrix axes name a parasite
        column_name = f'{parasite[0]}. {parasite.split(" ")[1]} proteins'
        table[column_name] = degree[parasite].fillna(0).astype(int).values
        columns.append(column_name)

    surface = summarise_localisations(edges, localisations)
    if surface is not None:
        table['DeepLoc'] = [surface['surface'].get(names[t], NO_LOCALISATION)
                            for t in degree.index]

    # the proteins reached by the most parasite proteins first
    table = table.assign(_reach=table[columns].sum(axis=1)).sort_values(
        ['_reach', 'Host protein'], ascending=[False, True], kind='stable').drop(columns='_reach')

    st.dataframe(table, width='stretch', hide_index=True)
    st.download_button('Download table', table.to_csv(index=False).encode('utf-8'),
                       file_name=f'shared_interactors_{first}_{second}.csv'.replace(' ', '_'),
                       mime='text/csv')


@st.cache_data(show_spinner=False)
def get_top_shared_proteins(df_pred, groups, group_order, niches, order_by,
                            annotations=None, localisations=None, top=40):
    '''The host proteins that the most parasites are predicted to interact with.'''
    edges = df_pred[['taxid1_label', 'source', 'target', 'target_name']].drop_duplicates()
    # one row per dot, keyed by display gene name so aliases do not stack
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
    '''
    called = dots.drop_duplicates('protein').set_index('protein')
    scores = {c: web_utils.DEEPLOC_SCORES[c] for c in web_utils.HOST_CLASSES
              if web_utils.DEEPLOC_SCORES[c] in called}
    labels = []
    for protein in proteins:
        row = called.loc[protein]
        crossed = [c for c, column in scores.items()
                   if row[column] > web_utils.DEEPLOC_CUTOFFS[c]]
        # a protein of one class is named by it, one with no localisation by what the legend
        # calls it
        labels.append(' + '.join(crossed) if len(crossed) > 1 else row['surface'])

    return labels


def add_parasite_strips(figure, dots, parasites, palette):
    '''The taxonomic group and the niche of each parasite, as two bands above its column.'''
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


def strip_width(parasites, proteins, height, room):
    '''
    Width in x units of the band beside the rows of the dot plot, at which its cells
    come out as wide as they are high.
    '''
    # the y axis is opened up past the top row to hold the two strips
    row = (height - DOT_CHROME) / (proteins + 2.5)
    # room = row * (parasites + 1.5 * width) / width, solved for the width
    return row * parasites / max(room - 1.5 * row, row)


def add_localisation_strip(figure, dots, proteins, width):
    '''Where DeepLoc puts each host protein, as a band beside its row.'''
    surface_of = dict(zip(dots['protein'], dots['surface']))
    localisations = [surface_of.get(p, NO_LOCALISATION) for p in proteins]
    shown = [c for c in LOCALISATION_ORDER if c in set(localisations)]
    scale, codes = flat_colour_scale(localisations,
                                     {c: LOCALISATION_COLORS[c] for c in shown})
    # the band is `width` x units across, an x unit being a parasite
    figure.add_trace(go.Heatmap(z=[[code] for code in codes], x0=-0.5 - width, dx=width,
                                y=list(range(len(proteins))),
                                text=[[label] for label in
                                      localisation_labels(dots, proteins)],
                                ygap=1, **scale))

    return shown


def add_dot_legend(figure, groups, localisations, niches, palette, width, height):
    '''
    Names the two strips the dot plot carries no legend for, and lays the three keys of
    the figure out under it as three legends of their own.
    '''
    # dots are sized by degree, so the legend is drawn from markers of one size instead of
    # the traces
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

    # measured from the foot of the figure, since the top margin plotly takes for the names
    # is not known here
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
    parasites in the order the other figures use so the taxonomic groups stay together,
    and the proteins ordered by how many parasites reach them.
    '''
    localised = 'surface' in dots.columns
    cells = {'x': {p: i for i, p in enumerate(parasites)},
             'y': {p: i for i, p in enumerate(proteins)}}
    dots = dots.assign(x=dots['parasite'].map(cells['x']),
                       y=dots['protein'].map(cells['y'])).dropna(subset=['x', 'y'])

    # written out: `parasites` and `degree` are two different counts, and the parasite is
    # named from the row rather than the x
    hover_columns = ['protein_full', 'parasites', 'degree', 'parasite']
    hover_lines = ['%{customdata[0]}', 'parasites reaching it: %{customdata[1]}',
                   'proteins of %{customdata[3]} reaching it: %{customdata[2]}']
    if localised:
        # the probabilities behind the band, two to a line
        scores = [web_utils.DEEPLOC_SCORES[c] for c in web_utils.HOST_CLASSES
                  if web_utils.DEEPLOC_SCORES[c] in dots]
        parts = [f'{DEEPLOC_LABELS[name]} %{{customdata[{len(hover_columns) + i}]:.2f}}'
                 for i, name in enumerate(scores)]
        hover_columns += scores
        hover_lines += [', '.join(parts[i:i + 2]) for i in range(0, len(parts), 2)]

    size = dot_size(parasites, proteins, column)
    # the area of the largest dot stands for the largest degree
    sizeref = max(1, dots['degree'].max()) / size ** 2
    figure = go.Figure()
    # one trace per group, so the groups are the legend
    for group in [g for g in list(palette) + [UNKNOWN_GROUP] if g in set(dots['group'])]:
        rows = dots[dots['group'] == group]
        figure.add_trace(go.Scatter(
            x=rows['x'], y=rows['y'], mode='markers', name=group,
            marker=dict(color=palette.get(group, UNKNOWN_COLOR), size=rows['degree'],
                        sizemode='area', sizeref=sizeref, sizemin=4, line=dict(width=0)),
            customdata=rows[hover_columns].to_numpy(),
            hovertemplate='<br>'.join(hover_lines) + '<extra></extra>'))

    groups_shown, niches = add_parasite_strips(figure, dots, parasites, palette)
    # the legends wrap inside what the row labels leave of the column
    height = max(420, DOT_ROW * len(proteins) + DOT_CHROME)
    room = plot_room(proteins, column)
    band = strip_width(len(parasites), len(proteins), height, room) if localised else 0
    localisations = add_localisation_strip(figure, dots, proteins, band) if localised else []
    legends = add_dot_legend(figure, groups_shown, localisations, niches, palette, room,
                             height)

    figure.update_layout(height=height, plot_bgcolor='white',
                         # the legends are placed against the foot of the figure, under the
                         # plot
                         margin=dict(l=0, r=0, t=10, b=legends),
                         xaxis_title=None, yaxis_title='host protein')
    # every parasite is named at the size a column has room for; the range runs back past
    # the band beside the rows
    figure.update_xaxes(range=[-0.5 - 1.5 * band, len(parasites) - 0.5],
                        side='top', tickmode='array', tickvals=list(range(len(parasites))),
                        ticktext=parasites, tickangle=-60, automargin=True, ticks='',
                        tickfont=dict(size=max(SMALLEST_LABEL, min(11, round(size / 1.1)))),
                        showgrid=True, gridcolor='#f0f0f0', zeroline=False)
    # reversed, so the most-shared protein is the top row
    figure.update_yaxes(range=[len(proteins) - 0.5, -3], tickmode='array',
                        tickvals=list(range(len(proteins))), ticktext=proteins, ticks='',
                        automargin=True, showgrid=True, gridcolor='#f0f0f0', zeroline=False)

    return figure


# the matrix keeps its cells square, so it needs the column width; drawn for
# DEFAULT_PAGE_WIDTH until the browser answers
column = web_utils.column_width(2)

st.caption('The parasites predicted against one host, compared with each other: which host '
           'interactors they share, which host proteins several of them reach, and the '
           'tissues and cell types in which their interactions can take place.')
st.markdown("---")

col1, col2, col3 = st.columns(3)

with col1:
    st.write('')

with col2:
    selected_host, selected_taxids = web_utils.host_selector(
        config, web_utils.load_predictions(data_dir),
        'Select a host to compare the parasites that infect it')
    if selected_host == web_utils.NO_HOST:
        st.text('Choose 1 host to explore the predicted host-parasite interactions')

with col3:
    st.write('')


if selected_host != web_utils.NO_HOST:
    # the heatmap and dot matrix read the same filtered interactions
    parasite_groups = {p['label']: p.get('group', UNKNOWN_GROUP)
                       for p in config['parasites'].values()}
    group_order = {g: i for i, g in enumerate(config.get('parasite_groups', {}))}
    niches = web_utils.get_niches(config)

    # one slider for the three figures below; the tissue dots at the foot keep every
    # prediction
    slider_column, order_column = st.columns([2, 1])
    with slider_column:
        score = st.slider('Confidence score', MIN_SCORE, MAX_SCORE, DEFAULT_SCORE,
                          help='Interactions predicted below this confidence are left out of '
                               'the three figures below. The tissue plot at the foot of the '
                               'page counts every prediction.')
    with order_column:
        # the order decides which blocks the figures show: contiguous values make a square
        # against the diagonal
        order_by = st.radio('Order the parasites by', [ORDER_BY_GROUP, ORDER_BY_NICHE],
                            horizontal=True,
                            help='Which annotation the two figures below put next to each '
                                 'other on their axes. The strips beside the axes show both '
                                 'either way; this is which of them comes out in blocks.')
    counted = get_tissue_expressed_predictions(data_dir, config, selected_taxids, score)
    shared_similarity = get_shared_interactor_similarity(counted, parasite_groups,
                                                        group_order, niches, order_by)

    per_tissue, per_cell_type = count_interactions_per_tissue(data_dir, config,
                                                              selected_taxids, score)
    ranked = per_tissue.groupby('Tissue')['interactions'].sum().sort_values(ascending=False,
                                                                           kind='stable')
    annotated = set(per_cell_type['Tissue'])
    # every tissue is offered, so a missing cell-type annotation is reported rather than
    # guessed at
    choices = list(ranked.index)

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
            # the figure carries the column width rather than being stretched, which keeps
            # the matrix square after full screen
            figure, click_layer = generate_shared_interactor_heatmap(
                *shared_similarity, config.get('parasite_groups', {}), column)
            # remounted after every dialog: Streamlit drops an identical selection and
            # plotly reads a second click as a deselection
            nonce = st.session_state.get(CELL_NONCE_KEY, 0)
            clicked = st.plotly_chart(
                figure, width='content', on_select='rerun', selection_mode='points',
                key=f'shared_cells_{selected_host}_{order_by}_{score}_{nonce}')
            # the layer the click came from is checked
            cell = next((point for point
                         in (clicked or {}).get('selection', {}).get('points', [])
                         if point.get('curve_number') == click_layer), None)
            if cell is not None:
                # the two parasites are the row and column the marker sits on
                parasites = list(shared_similarity[0].index)
                st.session_state[CELL_NONCE_KEY] = nonce + 1
                show_shared_interactors_dialog(
                    parasites[int(cell['y'])], parasites[int(cell['x'])], counted,
                    web_utils.load_protein_annotations(data_dir),
                    web_utils.load_deeploc_localisations(data_dir))
        else:
            st.text(f'Fewer than three parasites of {selected_host} share any host protein')

    with shared:
        st.subheader("Host interactors common to several parasites")
        # offered whatever the tissue leaves, so an emptying choice can be changed back
        shared_tissue = st.selectbox(
            'Tissue', [ALL_TISSUES] + choices, index=0, key='shared_tissue',
            help='Keep only the host proteins expressed in this tissue, for the parasites '
                 'known to infect it. Tissues the parasites infect, most interactions first.')
        # narrowed to one tissue where one is chosen; the heatmap keeps every tissue
        shared_predictions = (counted if shared_tissue == ALL_TISSUES else
                              get_tissue_expressed_predictions(data_dir, config,
                                                               selected_taxids, score,
                                                               shared_tissue))
        top_shared = get_top_shared_proteins(shared_predictions, parasite_groups, group_order,
                                             niches, order_by,
                                             web_utils.load_protein_annotations(data_dir),
                                             web_utils.load_deeploc_localisations(data_dir))
        where = '' if shared_tissue == ALL_TISSUES else f' in {shared_tissue}'
        if top_shared is not None:
            # the count of what the rows were taken from belongs to the caption
            *figure_arguments, shareable = top_shared
            shown = len(figure_arguments[1])
            # a truncated figure says what it is a top of
            selection = (f'The {shown} host proteins{where} reached by the most parasites, of '
                         f'the {shareable} reached by more than one.' if shareable > shown else
                         f'The {shown} host proteins{where} reached by more than one '
                         'parasite, most first.')
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
        else:
            st.info(f'No host protein{where} is reached by more than one parasite at this '
                    'confidence.')

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
                                  help='Tissues the parasites infect, most interactions first')
            if tissue in annotated:
                st.plotly_chart(generate_cell_type_bars(per_cell_type, tissue, parasite_groups,
                                                        config.get('parasite_groups', {})),
                                width='stretch')
            else:
                st.info(f'No single cell data available for {tissue} in {selected_host}, so '
                        'the interactions there cannot be split by cell type.')

st.markdown("---")


# Footer
with st.container():
    web_utils.footer()
