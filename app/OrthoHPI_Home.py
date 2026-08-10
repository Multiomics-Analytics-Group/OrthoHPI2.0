import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
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
hv.extension('bokeh')
from streamlit_bokeh import streamlit_bokeh

st.set_page_config(layout="wide", page_title="OrthoHPI 2.0", menu_items={})
st.session_state.data_dir = 'data'
st.session_state.config_file = 'config.yml'
style.load_css()

page = web_utils.show_pages_menu(index=0)
if page == "Predicted Host-parasite PPIs":
    st.switch_page("pages/1_Predicted_Host-Parasite_PPIs.py")
elif page == "Predicted PPI structures":
    st.switch_page('pages/2_Interaction_structures.py')
elif page == "About":
    st.switch_page('pages/3_About.py')
    

# Read dataset
config = utils.read_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()

# fallback for a parasite without a `group` in the config
UNKNOWN_GROUP = 'Unclassified'
UNKNOWN_COLOR = '#999999'
# colour of the chords at rest, before a parasite is hovered in the circos plot
REST_COLOR = '#bdbdbd'
# gradient the chords of the hovered parasite carry the number of shared host proteins
# on. Plasma reversed rather than a single-hue ramp: the steps of one hue were too close
# to read off the plot, and running through orange, magenta and purple to navy separates
# them by hue as well as by lightness. The pale yellow end is cut off, since a hairline
# in it is invisible against the white background, and the dark end is the many-proteins
# end so the pairs that share the most are the ones that stand out
CHORD_CMAP = ['#fb9f3a', '#ea7457', '#d04d73', '#ad2793', '#8104a7', '#4e02a2', '#0d0887']
# side of a cell of the similarity heatmap, in pixels, which is what sizes that figure
CELL = 22

#Initialize variables
df_select = None
net = None
selected_rows = []
selected_terms = []
enrichment_table = None
enrichment = None


@st.cache_data(show_spinner=False)
def load_predictions(data_dir):
    predictions = utils.read_parquet_file(input_file=f'{data_dir}/predictions.parquet')
    predictions['weight'] = predictions['weight'].astype(float)

    return predictions


@st.cache_data(show_spinner=False)
def get_host_predictions(data_dir, host_taxids):
    '''
    Predictions against one host group. host_taxids is a tuple of taxids (as strings)
    so grouped hosts -- rat and mouse under Rodent -- are pooled into a single view;
    the predictions keep their own species label and taxid.
    '''
    predictions = load_predictions(data_dir)

    return predictions[predictions['taxid2'].isin(host_taxids)]


@st.cache_data(show_spinner=False)
def count_interactions_per_tissue_cell_type(data_dir, config, host_taxids):
    '''
    Counts the predicted interactions of a host group per parasite, tissue and cell type,
    keeping only the tissues each parasite is known to infect (config['parasites']).
    Aggregating here keeps the merged predictions/tissues table out of the app's
    memory and makes sure each icicle sector is counted once.
    '''
    predictions = get_host_predictions(data_dir, host_taxids)
    tissues = utils.read_parquet_file(input_file=f'{data_dir}/tissues_cell_types.parquet')
    tpm_col = 'nTPM' if 'nTPM' in tissues.columns else 'pTPM'
    tissues = tissues.rename({'Gene': 'target'}, axis=1)[['target', 'Tissue', 'Cell type', tpm_col]]

    aux = pd.merge(predictions[['taxid1', 'taxid1_label', 'taxid2', 'target']], tissues, on='target', how='left')
    aux['taxid1'] = aux['taxid1'].astype(str)
    aux['Cell type'] = aux['Cell type'].fillna("Not available")

    mapped_tissues = config['tissues']
    infected_tissues = pd.DataFrame([(str(taxid), mapped_tissues[t].lower())
                                     for taxid, parasite in config['parasites'].items()
                                     for t in parasite['tissues']],
                                    columns=['taxid1', 'Tissue'])
    aux = pd.merge(aux, infected_tissues, on=['taxid1', 'Tissue'])

    counts = aux.groupby(['taxid1', 'taxid1_label', 'Tissue', 'Cell type'], observed=True).agg(
        edges_cell_type=('taxid2', 'count'), **{tpm_col: (tpm_col, 'mean')}).reset_index()
    counts['edges_tissue'] = counts.groupby(['taxid1', 'Tissue'])['edges_cell_type'].transform('sum')

    return counts


@st.cache_data(show_spinner=False)
def generate_tissue_cell_type_box(counts):
    tpm_col = 'nTPM' if 'nTPM' in counts.columns else 'pTPM'
    fig = px.icicle(counts, path=[px.Constant("Parasites"), 'taxid1_label', 'Tissue', 'Cell type'], values='edges_cell_type',
                  color='edges_cell_type', hover_data=['edges_tissue', 'edges_cell_type', 'taxid1', 'taxid1_label', tpm_col],
                  color_continuous_scale='Burgyl', height=900, maxdepth=-1)

    return fig

@st.cache_data(show_spinner=False)
def get_tissue_expressed_predictions(data_dir, config, host_taxids):
    '''
    Predictions restricted to the host proteins that are expressed in a tissue the
    parasite is known to infect (config['parasites'][taxid]['tissues']), which is the
    same restriction the tissue/cell type icicle applies. What is left are the
    interactions that could take place where the parasite actually is, rather than every
    interaction predicted from orthology. Parasites left without any interactor simply
    do not appear in what is built from this.
    '''
    predictions = get_host_predictions(data_dir, host_taxids)[['taxid1', 'taxid1_label',
                                                               'target', 'target_name']]
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

    return aux[['taxid1_label', 'target', 'target_name']].drop_duplicates()


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
def generate_circos_plot(data_dir, host_taxids, groups, palette, config):
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
    tissues it infects.
    '''
    predictions = get_tissue_expressed_predictions(data_dir, config, host_taxids)
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
    # ignored, which is what edge_line_color takes. On a linear ramp 236 of the 377 human
    # chords land in the first colour, since half the pairs share fewer than 5 proteins
    # and one shares 34, so the scale is logarithmic
    shading = ({'edge_color': dim('shared'), 'edge_cmap': CHORD_CMAP, 'logz': True} if graded
               else {'edge_line_color': CHORD_CMAP[-1]})
    chord = hv.Chord((links, hv.Dataset(nodes, 'index', ['name', 'group'])), vdims=['value', 'shared'])
    chord.opts(
        opts.Chord(width=500, height=700, labels='name',
                   node_color=dim('group').str(), cmap=palette, node_line_color='white',
                   **shading, edge_line_width=1.1,
                   edge_alpha=0.35, colorbar=graded,
                   colorbar_opts={'title': 'shared host proteins', 'title_text_font_style': 'normal',
                                  'title_text_font_size': '10px', 'width': 10, 'padding': 2},
                   edge_hover_line_alpha=1, edge_hover_line_width=1.8,
                   edge_selection_line_alpha=1, edge_selection_line_width=1.8,
                   edge_nonselection_line_alpha=0.1,
                   inspection_policy='nodes'))

    shown = set(nodes['group'])

    figure = separate_arcs(inset_chord_ends(gradient_on_hover(hv.render(chord))))

    return figure, [(g, c) for g, c in palette.items() if g in shown]

@st.cache_data(show_spinner=False)
def get_interactor_similarity(df_pred, groups, group_order):
    '''
    Jaccard similarity between the host proteins each pair of parasites is predicted to
    interact with. Jaccard rather than the count of shared proteins the circos draws: the
    parasites have between 12 and 200 predicted interactors, so a raw count mostly reports
    how many predictions a parasite has.

    The parasites are in the order of the circos and the dot matrix -- taxonomic group,
    then name -- so that a row is the same parasite in all three, and a clade is a block
    against the diagonal.
    '''
    targets = {g: set(df['target']) for g, df in df_pred.groupby('taxid1_label')}
    targets = {g: t for g, t in targets.items() if t}
    labels = sorted(targets, key=lambda p: (group_order.get(groups.get(p), len(group_order)), p))
    if len(labels) < 3:
        return None

    similarity = np.array([[len(targets[a] & targets[b]) / len(targets[a] | targets[b])
                            for b in labels] for a in labels])

    return (pd.DataFrame(similarity, index=labels, columns=labels),
            [groups.get(g, UNKNOWN_GROUP) for g in labels])


@st.cache_data(show_spinner=False)
def generate_similarity_heatmap(similarity, clades, palette):
    '''
    The similarity matrix, with a strip of the taxonomic group of each parasite down the
    side and along the bottom, so that the two axes are visibly the same list of parasites
    in the same order. The cells are held square (scaleanchor), which is the other half of
    reading the matrix as symmetric.

    Squaring the cells means one of the two axes has to give up whatever space the figure
    has beyond the square, which plotly takes off the range and turns into blank margins
    inside the plot. So the figure is sized to the aspect the square already needs -- the
    matrix at CELL pixels a side plus room for the labels around it -- and is drawn at that
    size instead of being stretched to the page. scaleanchor is then only making up the
    difference between the room the labels were given and the room they take.

    The diagonal carries its true value of 1 rather than being left blank. Blank renders as
    the white the colour scale starts at, so it could not be told from a pair sharing
    nothing -- and 238 of the 992 off-diagonal cells of the human matrix are exactly that.
    '''
    shown = [g for g in palette if g in set(clades)]
    steps = [[i / len(shown), palette[g]] for i, g in enumerate(shown)]
    steps += [[(i + 1) / len(shown), palette[g]] for i, g in enumerate(shown)]
    names = [f'{g[0]}. {g.split(" ")[1]}' for g in similarity.index]
    strip = dict(colorscale=sorted(steps), zmin=-0.5, zmax=len(shown) - 0.5, showscale=False,
                 hovertemplate='%{text}<extra></extra>')

    figure = make_subplots(rows=2, cols=2, column_widths=[0.03, 0.97], row_heights=[0.97, 0.03],
                           horizontal_spacing=0.01, vertical_spacing=0.012)
    figure.add_trace(go.Heatmap(z=[[shown.index(c)] for c in clades], y=names,
                                text=[[c] for c in clades], xgap=0, ygap=1, **strip), row=1, col=1)

    figure.add_trace(go.Heatmap(z=similarity.to_numpy(), x=names, y=names, colorscale='Blues',
                                zmin=0, zmax=1,
                                hovertemplate='%{y} and %{x}<br>Jaccard %{z:.2f}<extra></extra>',
                                colorbar=dict(title='Shared<br>interactors<br>(Jaccard)',
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
    figure.update_xaxes(tickangle=-60, showticklabels=True, row=2, col=2)
    figure.update_yaxes(showticklabels=False, row=2, col=2)
    figure.update_xaxes(visible=False, row=2, col=1)
    figure.update_yaxes(visible=False, row=2, col=1)
    figure.update_yaxes(autorange='reversed')

    # the widest parasite name, which is what the labels need down the left and, turned
    # through 60 degrees, under the bottom
    label = 6.5 * max(len(n) for n in names) + 12
    left, right, top, bottom = label, 130, 60, 0.87 * label + 25
    side = min(720, max(240, CELL * len(names)))
    figure.update_layout(width=side / 0.96 + left + right, height=side / 0.958 + top + bottom,
                         plot_bgcolor='white',
                         margin=dict(l=left, r=right, t=top, b=bottom),
                         legend=dict(orientation='h', yanchor='bottom', y=1.01, xanchor='left',
                                     x=0, itemclick=False, itemdoubleclick=False,
                                     font=dict(size=11)))

    return figure


@st.cache_data(show_spinner=False)
def get_top_shared_proteins(df_pred, groups, group_order, top=40):
    '''
    The host proteins that the most parasites are predicted to interact with. The circos
    and the heatmap both count how much two parasites have in common but neither says what
    they have in common, which is what this is for. Proteins only one parasite interacts
    with are left out: they are not shared by anything.
    '''
    pairs = df_pred[['taxid1_label', 'target', 'target_name']].drop_duplicates()
    counts = pairs.groupby('target_name')['taxid1_label'].nunique()
    counts = counts[counts > 1].sort_values(ascending=False, kind='stable')
    if counts.empty:
        return None

    proteins = list(counts.head(top).index)
    dots = pairs[pairs['target_name'].isin(proteins)].copy()
    dots['group'] = dots['taxid1_label'].map(lambda p: groups.get(p, UNKNOWN_GROUP))
    dots['parasites'] = dots['target_name'].map(counts)
    dots['parasite'] = dots['taxid1_label'].map(lambda p: f'{p[0]}. {p.split(" ")[1]}')
    order = sorted(dots['taxid1_label'].unique(),
                   key=lambda p: (group_order.get(groups.get(p), len(group_order)), p))

    return dots, proteins, [f'{p[0]}. {p.split(" ")[1]}' for p in order], int(counts.max())


@st.cache_data(show_spinner=False)
def generate_shared_protein_dots(dots, proteins, parasites, most, palette):
    '''
    A dot wherever a parasite is predicted to interact with one of the proteins, the
    parasites in the order of the circos so the taxonomic groups stay together, and the
    proteins ordered by how many parasites reach them.
    '''
    figure = px.scatter(dots, x='parasite', y='target_name', color='group',
                        # plotly express flips category_orders on a y axis, so `proteins`
                        # most-shared first puts the most-shared protein in the top row
                        color_discrete_map=palette, category_orders={
                            'parasite': parasites, 'target_name': proteins,
                            'group': [g for g in palette if g in set(dots['group'])]},
                        hover_data={'parasites': True, 'parasite': True, 'target_name': True,
                                    'group': False})
    figure.update_traces(marker=dict(size=8, line=dict(width=0)))
    figure.update_layout(height=max(420, 19 * len(proteins) + 240), plot_bgcolor='white',
                         margin=dict(l=0, r=0, t=10, b=10), legend_title_text='',
                         legend=dict(orientation='h', yanchor='bottom', y=1.01, x=0),
                         xaxis_title=None, yaxis_title=f'host protein (up to {most} parasites)')
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

def show_circos_plot(data_dir, host, host_taxids, config, caption, key):
    st.caption(caption)
    circos_plot, circos_legend = generate_circos_plot(
        data_dir, host_taxids,
        {p['label']: p.get('group', UNKNOWN_GROUP) for p in config['parasites'].values()},
        config.get('parasite_groups', {}), config)
    if circos_plot is None:
        st.text(f'No host proteins are shared by the parasites infecting {host} '
                'in the tissues they infect')
        return

    show_circos_legend(circos_legend)
    streamlit_bokeh(circos_plot, use_container_width=True, key=key)

def generate_boxplot_score_stats(df):
    fig = px.box(df.sort_values("taxid1"), x="taxid1_label", y="weight", color='taxid1', labels={"weight":"score", "taxid1_label": "parasites"})
    fig.update_traces(showlegend=False)

    return fig

def generate_barplot_stats(df):
    fig = px.bar(df.groupby(["taxid1_label"]).count().reset_index().sort_values("taxid2"), x="taxid1_label", y="weight", color='taxid1_label',  labels={"weight":"count", "taxid1_label": "parasites"})
    fig.update_traces(showlegend=False)
    return fig

@st.cache_data(show_spinner=False)
def generate_stats_plots(df):
    stats_figures = []
    fig = generate_barplot_stats(df)
    stats_figures.append((fig, "Number of Interactions",
                          "How many interactions with the host are predicted for each parasite."))

    fig = generate_boxplot_score_stats(df)
    stats_figures.append((fig, "Boxplot of Confidence scores",
                          "Distribution of the confidence score of each parasite's predicted "
                          "interactions: the higher the score, the better the evidence "
                          "transferred from the orthologous interaction."))

    return stats_figures




st.markdown("<h1 style='text-align: center; color: #023858;'>OrthoHPI 2.0</h1>", unsafe_allow_html=True)
st.markdown("<h3 style='text-align: center; color: #2b8cbe;'>Orthology Prediction of Host-Parasite PPI</h3>", unsafe_allow_html=True)

st.text(" ")
st.text(" ")
st.markdown("---")

col1, col2, col3 = st.columns(3)

with col1:
    st.write('')

with col2:
    # the host applies to every page (web_utils.host_selector), and hosts the config
    # groups together (rat + mouse) are one option
    selected_host, selected_taxids = web_utils.host_selector(
        config, load_predictions(data_dir), 'Select a host to visualize the predicted interactions')
    if selected_host == web_utils.NO_HOST:
        st.text('Choose 1 host to explore the predicted host-parasite interactions')

with col3:
    st.write('')


if selected_host != web_utils.NO_HOST:
    chart1, chart2 = st.columns(2)

    with chart1:
        st.subheader("Circos Plot of Common Host Interactors")
        caption = (f'Each arc is a parasite infecting {selected_host}, coloured by its taxonomic '
                   'group; a chord joins two parasites that are predicted to interact with the '
                   'same host proteins. Hover (or click) a parasite to pick out its chords, which '
                   'are then coloured by how many host proteins each pair shares, on the scale to '
                   'the right. Only the host proteins expressed in the tissues each parasite is '
                   'known to infect are counted, so what is left are the interactions that can '
                   'take place where the parasite is.')
        show_circos_plot(data_dir, selected_host, selected_taxids, config, caption, 'circos')

    # the heatmap and the dot matrix count the same interactions as the circos: the heatmap
    # reads where the circos cannot, since it normalises for how many predictions a parasite
    # has and has no chords to overlap, and the dot matrix names the host proteins that
    # neither of the other two ever shows
    parasite_groups = {p['label']: p.get('group', UNKNOWN_GROUP)
                       for p in config['parasites'].values()}
    group_order = {g: i for i, g in enumerate(config.get('parasite_groups', {}))}
    counted = get_tissue_expressed_predictions(data_dir, config, selected_taxids)
    similarity = get_interactor_similarity(counted, parasite_groups, group_order)
    top_shared = get_top_shared_proteins(counted, parasite_groups, group_order)

    stats_figs = generate_stats_plots(get_host_predictions(data_dir, selected_taxids))
    stats_cols = st.columns(len(stats_figs))
    i = 0
    for stats_fig, title, caption in stats_figs:
        with stats_cols[i]:
            st.subheader(title)
            st.caption(caption)
            st.plotly_chart(stats_fig, width='stretch')
        i += 1

    fig = generate_tissue_cell_type_box(count_interactions_per_tissue_cell_type(data_dir, config, selected_taxids))
    with chart2:
        st.subheader("Summary of Interactions per Tissue and Cell type")
        st.caption('Where the predicted interactions can take place: each parasite is broken down '
                   'into the tissues it is known to infect and the cell types of those tissues, '
                   'sized and coloured by the number of interactions. Click a block to zoom in.')
        st.plotly_chart(fig, width='stretch')

    if similarity is not None:
        st.subheader("Host Interactors Shared by Each Pair of Parasites")
        st.caption('How much of their predicted host interactors two parasites have in common, '
                   'as a share of everything either of them interacts with, so a parasite with '
                   'many predictions does not look similar to everything. The parasites are in '
                   'the order of the circos down both axes -- taxonomic group, then name, with '
                   'a strip of the group along each -- so a dark block against the diagonal is a '
                   'clade whose parasites converge on the same host proteins, and a dark cell '
                   'away from it is two unrelated parasites that do. The diagonal is each '
                   'parasite against itself and is 1 by definition.')
        # the figure carries its own size, since a square matrix stretched to the page is
        # a square with blank space either side of it rather than a wider square
        st.plotly_chart(generate_similarity_heatmap(*similarity, config.get('parasite_groups', {})),
                        width='content')

    if top_shared is not None:
        st.subheader("Which Host Proteins the Parasites Have in Common")
        st.caption('The host proteins reached by the most parasites, with a dot wherever a '
                   'parasite is predicted to interact with one of them. The parasites are in the '
                   'order of the circos, so a row of dots across a whole taxonomic group is a '
                   'protein that group converges on, and a gap is a parasite that does not reach '
                   'it. Proteins only one parasite interacts with are left out.')
        st.plotly_chart(generate_shared_protein_dots(*top_shared, config.get('parasite_groups', {})),
                        width='stretch')

st.markdown("---")


# Footer
with st.container():
    web_utils.footer()