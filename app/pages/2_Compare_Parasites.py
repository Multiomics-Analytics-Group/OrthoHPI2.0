import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import utils
import web_utils
import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import holoviews as hv
from plotly.subplots import make_subplots
from css import style
from holoviews import opts, dim
from bokeh.models import ColorBar, HoverTool, Text
hv.extension('bokeh')
from streamlit_bokeh import streamlit_bokeh

st.set_page_config(layout="wide", page_title="OrthoHPI 2.0", menu_items={})
style.load_css()
web_utils.show_header('Parasites of a host')

# Read dataset
config = utils.read_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()


# fallback for a parasite without a `group` in the config
UNKNOWN_GROUP = 'Unclassified'
UNKNOWN_COLOR = '#999999'
# colour of the chords at rest, before a parasite is hovered in the circos plot
REST_COLOR = '#bdbdbd'
# confidence the figures of shared interactors start at, and the range the slider spans.
# The same default and range as the network page, so a parasite shows the same
# interactions here as it does there
MIN_SCORE, MAX_SCORE, DEFAULT_SCORE = 0.4, 0.9, 0.4
# gradient the chords of the hovered parasite carry the number of shared host proteins
# on: light blue through to the dark blue of the page headings, which is also the dark
# end of the Blues the shared-interactor heatmap is drawn in. The palest steps of the ramp are
# left out -- a hairline in them is invisible against the white background -- and the
# dark end is the many-proteins end, so the pairs that share the most stand out
CHORD_CMAP = ['#a6bddb', '#74a9cf', '#3690c0', '#0570b0', '#045a8d', '#034368', '#023858']
# side of a cell of the shared-interactor heatmap, in pixels, which is what sizes that figure
CELL = 22
# how wide the two columns of the shared-interactors row come out, in pixels, on the window
# the app is being read in. Neither figure can ask its column how wide it is -- the matrix
# has to carry its own size to keep its cells square, and the circle has to carry its own
# height to stay round -- so both are sized from these. They are in the 1.2:1 ratio the
# columns are split in; raise them together for a wider window, lower them for a narrower
# one, and the figures follow.
MATRIX_COLUMN = 900
CIRCOS_COLUMN = 750
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
                   web_utils.NOT_SURFACE: 'square', NO_LOCALISATION: 'cross'}
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
    while the cell type counts keep only the cell types the protein is concentrated in
    (web_utils.keep_peak_cell_types, which says why). Host proteins the HPA gives no cell
    type -- every host but human, the single cell data being human only -- are counted in
    their tissue and left out of the cell type frame.
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
            pairs(web_utils.keep_peak_cell_types(aux), 'taxid1_label', 'Tissue', 'Cell type'))


@st.cache_data(show_spinner=False)
def generate_tissue_dots(per_tissue, groups, group_order, palette):
    '''
    A dot wherever a parasite is predicted to interact with a host protein expressed in a
    tissue it infects, sized by how many such interactions there are. The parasites are in
    the order of the circos and the two matrices above it, so a column is the same parasite
    everywhere on the page and a clade stays together.

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

    figure = px.scatter(dots, x='parasite', y='Tissue', color='group',
                        # plotly express flips category_orders on a y axis, so most-reached
                        # first puts the tissue the most parasites infect in the top row
                        color_discrete_map=palette, category_orders={
                            'parasite': parasites, 'Tissue': tissues,
                            'group': [g for g in palette if g in set(dots['group'])]},
                        size='interactions', size_max=18,
                        custom_data=['interactions'])
    figure.update_traces(marker=dict(sizemin=4, line=dict(width=0)),
                         hovertemplate='%{y}<br>%{x}<br>predicted interactions: '
                                       '%{customdata[0]}<extra></extra>')
    figure.update_layout(height=max(420, 19 * len(tissues) + 240), plot_bgcolor='white',
                         margin=dict(l=0, r=0, t=10, b=10), legend_title_text='',
                         legend=dict(orientation='h', yanchor='bottom', y=1.01, x=0),
                         xaxis_title=None, yaxis_title='tissue the parasite infects')
    figure.update_xaxes(tickangle=-60, showgrid=True, gridcolor='#f0f0f0')
    figure.update_yaxes(showgrid=True, gridcolor='#f0f0f0')

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

    A bar counts the interactions with host proteins concentrated in that cell type, which
    count_interactions_per_tissue defines; a protein abundant in several is counted in each,
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
    network page offers. It is applied here rather than per figure so the circos, the
    shared-interactor heatmap and the shared-protein dots keep counting the same interactions.

    One row is one predicted interaction, parasite protein (`source`) included: the circos
    and the heatmap only ever look at which host proteins are reached, but the dot matrix
    sizes its dots by how many parasite proteins reach each one.
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


@st.cache_data(show_spinner=False)
def get_common_interactors(df_pred, groups, group_order):
    '''
    Builds the circos nodes (parasites) and links (host proteins a pair of parasites
    shares). The parasites are ordered by taxonomic group and then by name so each
    clade occupies a contiguous arc of the circle, which is what makes the group
    colouring readable. The groups follow the order they are declared in
    config['parasite_groups'], so the circle and the legend read the same way.

    Only the pairs that share something are kept, and every link is given the same
    weight of 1: holoviews sizes the arc of a parasite from the weights of the chords
    that end on it, and the number of shared proteins is carried by `shared` for the
    colour to read instead.
    '''
    targets = {g: set(df['target']) for g, df in df_pred.groupby('taxid1_label')}
    labels = sorted(targets, key=lambda p: (group_order.get(groups.get(p), len(group_order)), p))

    nodes = [(i, g[0]+'. '+g.split(' ')[1], groups.get(g, UNKNOWN_GROUP)) for i, g in enumerate(labels)]
    links = [(i, j, 1, shared)
             for i, g1 in enumerate(labels) for j in range(i + 1, len(labels))
             if (shared := len(targets[g1].intersection(targets[labels[j]])))]

    return (pd.DataFrame(links, columns=['source', 'target', 'value', 'shared']),
            pd.DataFrame(nodes, columns=['index', 'name', 'group']))

@st.cache_resource(show_spinner=False)
def generate_circos_plot(data_dir, host_taxids, groups, palette, config, score=MIN_SCORE):
    '''
    The arcs are coloured by taxonomic group (config['parasite_groups']) rather than
    by species: a per-species palette would have to repeat itself for the 35 parasites
    infecting human, and the group colour also shows whether the shared host proteins
    follow the parasite phylogeny.

    The chords carry how many host proteins the two parasites share as a colour along
    CHORD_CMAP, all of them drawn at the same width. Sizing them instead is what a chord
    diagram normally does, but a pair sharing one protein and a pair sharing thirty are
    then a hairline and a band that swamps its neighbours, and with a chord for nearly
    every pair the thin end of the scale disappears. The colour is only shown for the
    parasite under the cursor (see gradient_on_hover) -- with all 377 chords of the human
    view coloured at once the circle is a wash that no scale can be read off. The hover
    styling has to set the alpha too: inheriting the resting alpha washes it out.

    The interactions are restricted to the host proteins the parasite could meet in the
    tissues it infects, and to those predicted at `score` confidence or better.
    '''
    predictions = get_tissue_expressed_predictions(data_dir, config, host_taxids, score)
    links, nodes = get_common_interactors(predictions, groups,
                                          {g: i for i, g in enumerate(palette)})
    if links.empty:
        return None, []

    palette = {**palette, **{g: UNKNOWN_COLOR for g in set(nodes['group']) if g not in palette}}
    # a scale needs something to range over: where every pair shares the same number of
    # proteins (Sus scrofa has a single chord) holoviews drops the mapping and falls back
    # to drawing the chords black, so those are given the top of the ramp and no colourbar
    graded = links['shared'].nunique() > 1
    # edge_color is the slot the colour dimension goes in and a literal colour there is
    # ignored, which is what edge_line_color takes. The scale is linear: it was logarithmic
    # for the counts EggNOG 5 gave, where half the pairs shared fewer than 5 proteins and a
    # linear ramp put 236 of the 377 human chords in the first colour. The EggNOG 6 counts
    # are far less skewed -- median 8 shared, max 77 -- so a linear ramp now spreads them
    shading = ({'edge_color': dim('shared'), 'edge_cmap': CHORD_CMAP} if graded
               else {'edge_line_color': CHORD_CMAP[-1]})
    chord = hv.Chord((links, hv.Dataset(nodes, 'index', ['name', 'group'])), vdims=['value', 'shared'])
    chord.opts(
        # the height is replaced by fit_circle_to_frame, which can only measure the labels
        # once they have been rendered; width is what streamlit stretches to the column
        opts.Chord(width=500, height=CIRCOS_COLUMN, labels='name',
                   node_color=dim('group').str(), cmap=palette, node_line_color='white',
                   **shading, edge_line_width=1.1,
                   edge_alpha=0.35, colorbar=graded,
                   # bokeh gives the bar the whole height of the plot unless it is told
                   # otherwise, which next to the circle is a stripe taller than the
                   # circle is wide for a scale of seven colours
                   colorbar_opts={'title': 'shared host proteins', 'title_text_font_style': 'normal',
                                  'width': 8, 'height': 160, 'padding': 2, 'label_standoff': 4},
                   # the title and the tick labels of the bar: holoviews writes its own
                   # sizes over anything colorbar_opts sets for them
                   fontsize={'clabel': '9pt', 'cticks': '9pt'},
                   edge_hover_line_alpha=1, edge_hover_line_width=1.8,
                   edge_selection_line_alpha=1, edge_selection_line_width=1.8,
                   edge_nonselection_line_alpha=0.1,
                   inspection_policy='nodes'))

    shown = set(nodes['group'])

    figure = fit_circle_to_frame(corner_colorbar(
        separate_arcs(inset_chord_ends(gradient_on_hover(hide_hover_tooltip(hv.render(chord)))))))

    return figure, [(g, c) for g, c in palette.items() if g in shown]


@st.cache_data(show_spinner=False)
def get_shared_interactor_counts(df_pred, groups, group_order):
    '''
    Number of host proteins each pair of parasites is predicted to interact with. The
    diagonal is fixed at the largest off-diagonal count as a visual boundary; every
    off-diagonal cell is the host proteins shared by that pair.

    The parasites are in the order of the circos and the dot matrix -- taxonomic group,
    then name -- so that a row is the same parasite in all three, and a clade is a block
    against the diagonal.
    '''
    targets = {g: set(df['target']) for g, df in df_pred.groupby('taxid1_label')}
    targets = {g: t for g, t in targets.items() if t}
    labels = sorted(targets, key=lambda p: (group_order.get(groups.get(p), len(group_order)), p))
    if len(labels) < 3:
        return None

    counts = np.array([[len(targets[a] & targets[b]) for b in labels] for a in labels])
    maximum_shared = np.triu(counts, k=1).max()
    np.fill_diagonal(counts, maximum_shared)

    return (pd.DataFrame(counts, index=labels, columns=labels),
            [groups.get(g, UNKNOWN_GROUP) for g in labels])


@st.cache_data(show_spinner=False)
def generate_shared_interactor_heatmap(counts, clades, palette, column=MATRIX_COLUMN):
    '''
    The shared-interactor count matrix, with a strip of the taxonomic group of each
    parasite down the side and along the bottom, so that the two axes are visibly the
    same list of parasites in the same order. The cells are held square (scaleanchor),
    which is the other half of reading the matrix as symmetric.

    Squaring the cells means one of the two axes has to give up whatever space the figure
    has beyond the square, which plotly takes off the range and turns into blank margins
    inside the plot. So the figure is sized to the aspect the square already needs -- the
    matrix at CELL pixels a side plus room for the labels around it -- and is drawn at that
    size instead of being stretched to the page. scaleanchor is then only making up the
    difference between the room the labels were given and the room they take.

    The diagonal carries the maximum shared-interactor count rather than being left
    blank. Blank renders as the white the colour scale starts at, so it could not be
    told from a pair sharing nothing.

    :param column: pixels the column holding this figure is expected to be. The matrix takes
                   what is left of it once the labels and the colour bar have taken theirs,
                   so the figure comes out the width of the column instead of sitting in the
                   middle of it -- and never wider than a cell of CELL pixels needs, since a
                   matrix of four parasites drawn across a whole column is four huge squares.
    '''
    shown = [g for g in palette if g in set(clades)]
    steps = [[i / len(shown), palette[g]] for i, g in enumerate(shown)]
    steps += [[(i + 1) / len(shown), palette[g]] for i, g in enumerate(shown)]
    names = [f'{g[0]}. {g.split(" ")[1]}' for g in counts.index]
    strip = dict(colorscale=sorted(steps), zmin=-0.5, zmax=len(shown) - 0.5, showscale=False,
                 hovertemplate='%{text}<extra></extra>')

    figure = make_subplots(rows=2, cols=2, column_widths=[0.03, 0.97], row_heights=[0.97, 0.03],
                           horizontal_spacing=0.01, vertical_spacing=0.012)
    figure.add_trace(go.Heatmap(z=[[shown.index(c)] for c in clades], y=names,
                                text=[[c] for c in clades], xgap=0, ygap=1, **strip), row=1, col=1)

    figure.add_trace(go.Heatmap(z=counts.to_numpy(), x=names, y=names, colorscale='Blues',
                                zmin=0, zmax=counts.to_numpy().max(),
                                hovertemplate='%{y} and %{x}<br>Shared interactors %{z:.0f}<extra></extra>',
                                colorbar=dict(title='Shared<br>interactors',
                                              thickness=12, len=0.6, y=1, yanchor='top')),
                     row=1, col=2)

    figure.add_trace(go.Heatmap(z=[[shown.index(c) for c in clades]], x=names,
                                text=[list(clades)], xgap=1, ygap=0, **strip), row=2, col=2)

    # the strips are heatmaps and cannot carry a legend of their own, so the groups are named
    # by empty traces whose only purpose is their legend entry
    for group in shown:
        figure.add_trace(go.Scatter(x=[None], y=[None], mode='markers', name=group,
                                    marker=dict(size=10, symbol='square', color=palette[group]),
                                    hoverinfo='skip', showlegend=True), row=1, col=2)

    # square cells, so the matrix reads as the symmetric thing it is whatever width the
    # browser gives it. The two strips follow the axes of the matrix rather than the other
    # way around: keeping `matches` off the axis that carries the constraint means that
    # when plotly pads a range to square the cells, the strips inherit the padding and stay
    # lined up with the rows and columns they label
    figure.update_yaxes(scaleanchor='x2', scaleratio=1, row=1, col=2)
    figure.update_yaxes(matches='y2', row=1, col=1)
    figure.update_xaxes(matches='x2', row=2, col=2)

    figure.update_xaxes(showticklabels=False, row=1, col=1)
    figure.update_xaxes(showticklabels=False, row=1, col=2)
    figure.update_yaxes(showticklabels=False, row=1, col=2)
    # every parasite is named on both axes, however small the cells are drawn: plotly thins
    # tick labels that no longer fit, and a heatmap with every third row labelled cannot be
    # read at all
    figure.update_xaxes(tickangle=-60, showticklabels=True, tickmode='linear', dtick=1,
                        tickfont=dict(size=10), row=2, col=2)
    figure.update_yaxes(showticklabels=False, row=2, col=2)
    figure.update_yaxes(tickmode='linear', dtick=1, tickfont=dict(size=10), row=1, col=1)
    figure.update_xaxes(visible=False, row=2, col=1)
    figure.update_yaxes(visible=False, row=2, col=1)
    figure.update_yaxes(autorange='reversed')

    # the widest parasite name, which is what the labels need down the left and, turned
    # through 60 degrees, under the bottom
    label = 6.5 * max(len(n) for n in names) + 12
    left, right, top, bottom = label, 105, 60, 0.87 * label + 25
    # 0.96 is the width the matrix is given of what is left, the rest going to the group
    # strip beside it and the gap between the two, so the figure comes out at `column`
    side = max(240, min((column - left - right) * 0.96, CELL * len(names)))
    figure.update_layout(width=side / 0.96 + left + right, height=side / 0.958 + top + bottom,
                         plot_bgcolor='white',
                         margin=dict(l=left, r=right, t=top, b=bottom),
                         legend=dict(orientation='h', yanchor='bottom', y=1.01, xanchor='left',
                                     x=0, itemclick=False, itemdoubleclick=False,
                                     font=dict(size=11)))

    return figure


# longest descriptive protein name written beside a gene symbol on the dot matrix. The
# label is an axis tick and every character of it is taken off the width of the plot, so
# the description is cut where it stops earning the space -- the whole of it is in the
# hover, which is where a name that long is read anyway.
DESCRIPTION_WIDTH = 38


def label_proteins(df_pred, annotations):
    '''
    Names each host protein by its gene symbol and its descriptive protein name, since the
    symbol on its own identifies the protein only for someone who already knows it.

    The descriptions are keyed by STRING id and the matrix is keyed by symbol -- a host
    group covering two species carries the same gene under an id of each -- so the first
    description found for a symbol is the one used. Proteins UniProt has nothing but
    "Uncharacterized protein" for keep their symbol alone: that description names nothing
    and would be repeated down the axis.

    :param df_pred: tissue-expressed predictions of the host group
    :param dict annotations: STRING id --> descriptive protein name
    :return: {gene symbol: axis label}
    '''
    described = {}
    for protein, name in df_pred[['target', 'target_name']].drop_duplicates().values:
        description = (annotations or {}).get(protein, '')
        if description and not description.lower().startswith('uncharacterized'):
            described.setdefault(name, description.rstrip('.'))

    labels = {}
    for name in df_pred['target_name'].unique():
        description = described.get(name)
        if description and len(description) > DESCRIPTION_WIDTH:
            description = description[:DESCRIPTION_WIDTH - 1].rstrip() + '…'
        labels[name] = f'{name} · {description}' if description else str(name)

    return labels


def summarise_localisations(df_pred, localisations):
    '''
    The DeepLoc call of each host protein of the matrix, keyed by the gene symbol the
    rows are drawn under rather than by the STRING id DeepLoc wrote it against.

    A host protein reaches this page only because DeepLoc called it surface-exposed, but
    the two ways of being surface-exposed are not the same interaction: a cell membrane
    protein is met on the surface of the host cell, an extracellular one is met in the
    fluid around it. Which of the two a protein was kept for is web_utils.classify_surface,
    which the home page reads the parasite proteins with.

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
def get_top_shared_proteins(df_pred, groups, group_order, annotations=None,
                            localisations=None, top=40):
    '''
    The host proteins that the most parasites are predicted to interact with. The circos
    and the heatmap both count how much two parasites have in common but neither says what
    they have in common, which is what this is for. Proteins only one parasite interacts
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
    labels = label_proteins(df_pred, annotations)
    dots['protein'] = dots['target_name'].map(labels)
    surface = summarise_localisations(df_pred, localisations)
    if surface is not None:
        dots = dots.join(surface, on='target_name')
        dots['surface'] = dots['surface'].fillna(NO_LOCALISATION)
        dots['localizations'] = dots['localizations'].fillna('')
    order = sorted(dots['taxid1_label'].unique(),
                   key=lambda p: (group_order.get(groups.get(p), len(group_order)), p))

    return (dots, [labels[p] for p in proteins],
            [f'{p[0]}. {p.split(" ")[1]}' for p in order], int(counts.max()))


def split_dot_legend(figure, groups, surfaces, palette):
    '''
    Splits the legend of the shared-interactors matrix back into its two channels.

    Plotly express draws a trace per combination of the two channels it is given and names
    it after both, so a matrix coloured by taxonomic group and shaped by surface class
    comes out with a legend of "Nematoda, Cell membrane" entries, one per combination that
    occurs. Every drawn trace is taken out of the legend and the legend is built from
    entries carrying no data instead: one per taxonomic group, in its colour and always as
    a circle so the shape of whichever combination came first says nothing, and one per
    surface class, in grey since the shape means the same whatever the colour it is drawn
    in. The group entries keep the legendgroup of the traces they name, so clicking one
    still hides that parasite group.

    :param figure: the matrix figure, modified in place
    :param list groups: taxonomic groups that occur, in the order of the legend
    :param list surfaces: surface classes that occur, in the order of the legend
    :param dict palette: {taxonomic group: colour}
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


@st.cache_data(show_spinner=False)
def generate_shared_protein_dots(dots, proteins, parasites, most, palette):
    '''
    A dot wherever a parasite is predicted to interact with one of the proteins, the
    parasites in the order of the circos so the taxonomic groups stay together, and the
    proteins ordered by how many parasites reach them. The dot is sized by how many
    proteins of that parasite reach that host protein.

    The area of the dot is what carries the degree (plotly's default), since that is the
    channel size is read on, and sizemin keeps the single-protein dots -- the largest group
    of them -- from collapsing to a speck next to a degree of eighty-six.

    Where the dots carry a DeepLoc call the shape of the marker is the surface class of
    the host protein, which is constant down a row: a circle is met on the cell membrane,
    a diamond in the fluid around the cell. The hover carries the two probabilities behind
    that and everywhere else DeepLoc puts the protein.
    '''
    localised = 'surface' in dots.columns
    orders = {'parasite': parasites, 'protein': proteins,
              'group': [g for g in palette if g in set(dots['group'])]}
    # the hover is written out rather than left to plotly express, which prints the raw
    # column name of whatever it is given. The protein and the parasite are the two axes
    # already, and `parasites` and `degree` are two different counts of two different
    # things, which as bare numbers under their column names they do not say
    hover_columns = ['parasites', 'degree']
    hover_lines = ['%{y}', 'parasites reaching it: %{customdata[0]}',
                   'proteins of %{x} reaching it: %{customdata[1]}']
    if localised:
        orders['surface'] = [s for s in SURFACE_SYMBOLS if s in set(dots['surface'])]
        hover_columns += ['cell_membrane', 'extracellular']
        # the class is the shape of the dot, so the hover carries the two probabilities
        # behind it rather than naming it a second time
        hover_lines.append('P(cell membrane) %{customdata[2]:.2f}, '
                           'P(extracellular) %{customdata[3]:.2f}')

    figure = px.scatter(dots, x='parasite', y='protein', color='group',
                        # plotly express flips category_orders on a y axis, so `proteins`
                        # most-shared first puts the most-shared protein in the top row
                        color_discrete_map=palette, category_orders=orders,
                        symbol='surface' if localised else None,
                        symbol_map=SURFACE_SYMBOLS if localised else {},
                        size='degree', size_max=15, labels=DEEPLOC_LABELS,
                        custom_data=hover_columns)
    figure.update_traces(marker=dict(sizemin=4, line=dict(width=0)),
                         hovertemplate='<br>'.join(hover_lines) + '<extra></extra>')
    if localised:
        split_dot_legend(figure, orders['group'], orders['surface'], palette)
    figure.update_layout(height=max(420, 19 * len(proteins) + 240), plot_bgcolor='white',
                         margin=dict(l=0, r=0, t=10, b=10), legend_title_text='',
                         legend=dict(orientation='h', yanchor='bottom', y=1.01, x=0),
                         xaxis_title=None, yaxis_title=f'host protein (up to {most} parasites)')
    # the labels are as long as a protein name, so plotly is left to take the room they
    # need off the plot rather than drawing them over it or cutting them at the margin
    figure.update_yaxes(automargin=True)
    figure.update_xaxes(tickangle=-60, showgrid=True, gridcolor='#f0f0f0')
    figure.update_yaxes(showgrid=True, gridcolor='#f0f0f0')

    return figure


def arc_angles(data):
    '''Start, end and middle angle of every node arc of a rendered chord.'''
    spans = []
    for xs, ys in zip(data['arc_xs'], data['arc_ys']):
        angles = np.unwrap(np.arctan2(ys, xs))
        spans.append((angles[0], angles[-1], (angles[0] + angles[-1]) / 2))

    return spans

def inset_chord_ends(figure, keep=0.86):
    '''
    Holoviews spreads the ends of a parasite's chords over its arc with a linspace that
    includes both ends, so the outermost chord of a parasite starts exactly on the boundary
    with its neighbour. Among the hundreds of chords of the human view that is invisible,
    but where a parasite carries two or three -- the other hosts -- a chord appears to run
    from one boundary to another instead of landing on a parasite at all. Pull every chord
    end towards the middle of the arc it belongs to and redraw the spline holoviews draws:
    a cubic bezier whose control points sit at half the radius of its ends.
    '''
    for renderer in figure.renderers:
        if not hasattr(renderer, 'node_renderer'):
            continue

        spans = arc_angles(renderer.node_renderer.data_source.data)
        edges = renderer.edge_renderer.data_source.data
        all_xs, all_ys = [], []
        for xs, ys, source, target in zip(edges['xs'], edges['ys'], edges['start'], edges['end']):
            xs, ys = np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)
            breaks = np.where(np.isnan(xs))[0]
            path_xs, path_ys = [], []
            for strand in np.split(np.arange(len(xs)), breaks):
                strand = strand[~np.isnan(xs[strand])]
                if len(strand) < 2:
                    continue
                ends = []
                for point, node in ((strand[0], source), (strand[-1], target)):
                    start, end, middle = spans[node]
                    angle = np.arctan2(ys[point], xs[point])
                    angle += 2 * np.pi * np.round((middle - angle) / (2 * np.pi))
                    angle = middle + (angle - middle) * keep
                    ends.append(np.array([np.cos(angle), np.sin(angle)]))

                steps = np.linspace(0, 1, len(strand))[:, None]
                curve = ((1 - steps) ** 3 * ends[0] + 3 * (1 - steps) ** 2 * steps * (ends[0] / 2)
                         + 3 * (1 - steps) * steps ** 2 * (ends[1] / 2) + steps ** 3 * ends[1])
                if path_xs:
                    path_xs.append(np.nan)
                    path_ys.append(np.nan)
                path_xs.extend(curve[:, 0])
                path_ys.extend(curve[:, 1])

            all_xs.append(np.array(path_xs))
            all_ys.append(np.array(path_ys))
        edges['xs'], edges['ys'] = all_xs, all_ys

    return figure

def gradient_on_hover(figure):
    '''
    Holds the chords at a neutral grey until a parasite is hovered or clicked, and gives
    the colour scale to the chords of that parasite only. Colouring all of them at rest
    puts every pair of the circle on screen at once, which is a wash of crossing lines
    where no chord can be followed to its ends, let alone have its colour compared with
    another; a dozen chords lit at a time can be read.

    Bokeh keeps a separate glyph for each state of a renderer, so this swaps what
    holoviews put on them: the colour mapping (a Field over the shared counts, or the
    flat colour where there is nothing to grade) moves to the hover and selection glyphs,
    and the resting and non-selected ones go grey.
    '''
    for renderer in figure.renderers:
        edges = getattr(renderer, 'edge_renderer', None)
        if edges is None:
            continue

        shading = edges.glyph.line_color
        for glyph in (edges.glyph, edges.nonselection_glyph, edges.muted_glyph):
            if glyph is not None:
                glyph.line_color = REST_COLOR
        for glyph in (edges.hover_glyph, edges.selection_glyph):
            if glyph is not None:
                glyph.line_color = shading

    return figure

def hide_hover_tooltip(figure):
    '''
    Drops the box the hover tool pops up over a parasite's arc and label. The tool itself
    stays -- it is what lights up the chords of the parasite under the cursor (see
    gradient_on_hover) -- and the box only repeated the name already written beside the
    arc, while covering the chords the hover had just brought out.
    '''
    for tool in figure.select(HoverTool):
        tool.tooltips = None

    return figure

def corner_colorbar(figure):
    '''
    Moves the scale of shared host proteins into the top right corner of the plot. Bokeh
    puts a colourbar in a panel beside the frame, which reserves a column of the figure
    for a bar 8 pixels wide and squeezes the circle for the rest of its height; the
    corners inside the frame are empty, since a circle drawn in a taller-than-wide figure
    leaves them, so the bar costs nothing there.
    '''
    for bar in figure.select(ColorBar):
        figure.right = [panel for panel in figure.right if panel is not bar]
        figure.add_layout(bar)
        bar.location = 'top_right'
        bar.orientation = 'vertical'
        # the panel painted it on the white of the figure; inside the frame the default
        # white would be a card sitting over the plot
        bar.background_fill_alpha = 0

    return figure

def fit_circle_to_frame(figure, column=CIRCOS_COLUMN, margin=1.05, floor=1.12):
    '''
    Fits the figure to what is actually drawn in it, top and bottom, which is not the same
    as the circle: the labels are written outwards along their own radius, so the label of
    the parasite at the top of the circle stands straight up and is the highest thing in the
    figure, while the label at three o'clock lies flat and reaches furthest sideways.

    Holoviews pads the ranges to 1.4 in every direction whatever the labels are -- 40% of a
    circle whose arcs stop at 1.05 -- so how much of that room is used depends on where the
    long names land. In the human view the labels reach 1.28 of it and the rest is white; a
    view of two parasites puts one label straight up and needs 1.42, more than it is given,
    and the name is cut off. It also leaves a title on the figure with nothing in it, which
    is an empty row about 23px tall above everything.

    So the labels are measured, the vertical range is set to what they reach, the height is
    set to draw that range at the same scale as the width, and the empty title is dropped.
    match_aspect keeps the circle round whatever width the column really turns out to be:
    the ranges are a starting point once it is on, and if the window is narrower than
    CIRCOS_COLUMN the circle is drawn smaller rather than squashed.

    :param column: pixels the figure is expected to be drawn at, which is what turns the
                   length of a label in pixels into the length of a label on these axes
    :param margin: what to leave beyond the last glyph, as a fraction
    :param floor: never crop closer to the circle than this, whatever the labels measure
    '''
    figure.title = None
    span = figure.x_range.end - figure.x_range.start
    per_unit = column / span
    reach = floor

    labels = next((r for r in figure.renderers
                   if isinstance(getattr(r, 'glyph', None), Text)), None)
    if labels is not None:
        data = labels.data_source.data
        # 0.52 em is about the average width of a glyph of the sans-serif bokeh sets text in
        size = float(str(labels.glyph.text_font_size).removesuffix('pt')) * 96 / 72
        radius = np.hypot(np.asarray(data['x'], dtype=float), np.asarray(data['y'], dtype=float))
        length = np.array([len(str(t)) * 0.52 * size for t in data['text']]) / per_unit
        tip = (radius + length) * np.sin(np.asarray(data['angle'], dtype=float))
        reach = max(floor, margin * np.abs(tip).max())

    figure.y_range.start, figure.y_range.end = -reach, reach
    figure.match_aspect = True
    figure.height = int(round(2 * reach * per_unit)) + 20   # the top and bottom borders

    return figure

def separate_arcs(figure):
    '''
    Holoviews draws each node arc from where the previous one ends, with no gap, so now
    that neighbouring parasites of the same taxonomic group share a colour they merge into
    one band and it is impossible to tell which of two neighbours a chord lands on. Mark
    every boundary with a white tick instead of shortening the arcs: a chord attaches
    anywhere along the arc of its parasite, so an arc that no longer covers its own angles
    leaves the outermost chords starting in mid air.
    '''
    for renderer in figure.renderers:
        data = getattr(renderer, 'data_source', None) and renderer.data_source.data
        if not data or 'arc_xs' not in data:
            continue

        angles = [np.arctan2(ys[-1], xs[-1]) for xs, ys in zip(data['arc_xs'], data['arc_ys'])]
        inner, outer = 0.972, 1.028
        figure.segment(x0=[inner * np.cos(a) for a in angles], y0=[inner * np.sin(a) for a in angles],
                       x1=[outer * np.cos(a) for a in angles], y1=[outer * np.sin(a) for a in angles],
                       line_color='white', line_width=2)
        break

    return figure

def show_circos_legend(legend):
    swatches = ''.join(
        f"<span style='display:inline-flex;align-items:center;margin:0 12px 4px 0;white-space:nowrap;'>"
        f"<span style='width:12px;height:12px;border-radius:2px;background:{color};"
        f"display:inline-block;margin-right:5px;'></span>{group}</span>"
        for group, color in legend)
    st.markdown(f"<div style='font-size:0.85em;color:#333333;'>{swatches}</div>", unsafe_allow_html=True)

def show_circos_plot(data_dir, host, host_taxids, config, caption, key, score):
    st.caption(caption)
    circos_plot, circos_legend = generate_circos_plot(
        data_dir, host_taxids,
        {p['label']: p.get('group', UNKNOWN_GROUP) for p in config['parasites'].values()},
        config.get('parasite_groups', {}), config, score)
    if circos_plot is None:
        st.text(f'No host proteins are shared by the parasites infecting {host} '
                f'in the tissues they infect, at confidence {score:g} or better')
        return

    show_circos_legend(circos_legend)
    streamlit_bokeh(circos_plot, use_container_width=True, key=key)


st.caption('The parasites predicted against one host, compared with each other: how '
           'similar their predicted host interactors are, which host proteins several of '
           'them reach, and the tissues and cell types in which their interactions can '
           'take place.')
st.markdown("---")

col1, col2, col3 = st.columns(3)

with col1:
    st.write('')

with col2:
    # the host applies to every page (web_utils.host_selector), and hosts the config
    # groups together (rat + mouse) are one option
    selected_host, selected_taxids = web_utils.host_selector(
        config, web_utils.load_predictions(data_dir),
        'Select a host to compare the parasites that infect it')
    if selected_host == web_utils.NO_HOST:
        st.text('Choose 1 host to explore the predicted host-parasite interactions')

with col3:
    st.write('')


if selected_host != web_utils.NO_HOST:
    # the heatmap and the dot matrix count the same interactions as the circos: the heatmap
    # reads where the circos cannot, since it normalises for how many predictions a parasite
    # has and has no chords to overlap, and the dot matrix names the host proteins that
    # neither of the other two ever shows
    parasite_groups = {p['label']: p.get('group', UNKNOWN_GROUP)
                       for p in config['parasites'].values()}
    group_order = {g: i for i, g in enumerate(config.get('parasite_groups', {}))}

    # one slider for the three figures below it, which are the same interactions read three
    # ways -- filtering them apart would put different numbers in figures whose captions
    # say they match. The tissue dots at the foot of the page are a different count and
    # keep every prediction.
    # It sits in the middle of three columns, as the host selector above it does: left to
    # itself a slider takes the whole width of the page, which is a metre of track for a
    # range of half a point
    with st.columns(3)[1]:
        score = st.slider('Confidence score', MIN_SCORE, MAX_SCORE, DEFAULT_SCORE,
                          help='Interactions predicted below this confidence are left out of '
                               'the three figures below. The tissue plot at the foot of the '
                               'page counts every prediction.')
    counted = get_tissue_expressed_predictions(data_dir, config, selected_taxids, score)
    shared_counts = get_shared_interactor_counts(counted, parasite_groups, group_order)
    top_shared = get_top_shared_proteins(counted, parasite_groups, group_order,
                                         web_utils.load_protein_annotations(data_dir),
                                         web_utils.load_deeploc_localisations(data_dir))

    # the two figures of how much the parasites share sit side by side, since they are the
    # same numbers read two ways: the matrix gives every pair a cell, the circle gives the
    # pairs that share anything a chord. The matrix is given the wider column -- it carries
    # a label per parasite on both of its axes, where the circle carries one per arc.
    matrix, circle = st.columns([1.2, 1])

    with matrix:
        st.subheader("Host interactors shared by each pair of parasites")
        st.caption('Number of host interactors shared by each pair of parasites. A strip of '
                   'the taxonomic group runs along each axis. The diagonal is fixed to the '
                   'largest shared-interactor count.')
        if shared_counts is not None:
            # the figure carries its own size, since a square matrix stretched to the page is
            # a square with blank space either side of it rather than a wider square
            st.plotly_chart(generate_shared_interactor_heatmap(
                *shared_counts, config.get('parasite_groups', {})),
                            width='content')
        else:
            st.text(f'Fewer than three parasites of {selected_host} share any host protein, '
                    'which is not a matrix worth drawing')

    with circle:
        st.subheader("Circos plot of common host interactors")
        caption = (f'Each arc is a parasite infecting {selected_host}, coloured by its taxonomic '
                   'group; a chord joins two parasites that are predicted to interact with the '
                   'same host proteins. Hover (or click) a parasite to pick out its chords, which '
                   'are then coloured by how many host proteins each pair shares, on the scale '
                   'in the corner.')
        show_circos_plot(data_dir, selected_host, selected_taxids, config, caption, 'circos', score)

    if top_shared is not None:
        st.subheader("Host interactors common to several parasites")
        st.caption('Host proteins reached by the most parasites, with a dot wherever a '
                   'parasite is predicted to interact with one, sized by the number of that '
                   "parasite's proteins reaching it. Proteins reached by a single parasite "
                   'are omitted. Dot shape gives the DeepLoc 2 localization of the host '
                   'protein, circles cell membrane and diamonds extracellular; hover for the '
                   'underlying probabilities.')
        st.plotly_chart(generate_shared_protein_dots(*top_shared, config.get('parasite_groups', {})),
                        width='stretch')

    st.subheader("Tissues in which the predicted interactions can take place")
    st.caption('Predicted interactions per parasite and tissue, restricted to the tissues '
               'each parasite is known to infect and sized by the number of interactions with '
               'proteins expressed there. Tissues are ordered by the number of parasites '
               'infecting them. An interaction is counted once per tissue, irrespective of '
               'the number of cell types the host protein is expressed in.')
    per_tissue, per_cell_type = count_interactions_per_tissue(data_dir, config,
                                                              selected_taxids, score)
    if per_tissue.empty:
        st.info('No predicted interaction is left at this confidence in a tissue the '
                'parasites are known to infect.')
    else:
        st.plotly_chart(generate_tissue_dots(per_tissue, parasite_groups, group_order,
                                             config.get('parasite_groups', {})),
                        width='stretch')

    # cell types are only worth drawing one tissue at a time, and only for a host they are
    # annotated for: the HPA single cell data is human, so every other host has no cell type
    # at all and nothing to draw
    ranked = per_tissue.groupby('Tissue')['interactions'].sum().sort_values(ascending=False,
                                                                           kind='stable')
    annotated = set(per_cell_type['Tissue'])
    choices = [t for t in ranked.index if t in annotated]
    if choices:
        st.subheader("Cell types of a tissue")
        st.caption('Predicted interactions per cell type of the selected tissue, stacked by '
                   'taxonomic group. A host protein counts towards a cell type where it '
                   'reaches at least half the expression it has anywhere in the tissue, so '
                   'the bars are the cell types the interaction partners are concentrated '
                   'in. A protein abundant in several counts in each, and the bars are '
                   'therefore not a partition of the tissue.')
        tissue = st.selectbox('Tissue', choices, index=0,
                              help='Tissues with cell type annotation, most interactions first')
        st.plotly_chart(generate_cell_type_bars(per_cell_type, tissue, parasite_groups,
                                                config.get('parasite_groups', {})),
                        width='stretch')

st.markdown("---")



# Footer
with st.container():
    web_utils.footer()
