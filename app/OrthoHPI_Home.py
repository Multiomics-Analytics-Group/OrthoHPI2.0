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
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform
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
# colour of the chords of the parasite under the cursor in the circos plot
HOVER_COLOR = '#00A000'
# circos view restricted to the host proteins expressed where the parasite is
IN_TISSUE_VIEW = 'In the infected tissues'

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
    predictions = get_host_predictions(data_dir, host_taxids)[['taxid1', 'taxid1_label', 'target']]
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

    return aux[['taxid1_label', 'target']].drop_duplicates()


@st.cache_data(show_spinner=False)
def get_common_interactors(df_pred, groups, group_order):
    '''
    Builds the circos nodes (parasites) and links (host proteins a pair of parasites
    shares). The parasites are ordered by taxonomic group and then by name so each
    clade occupies a contiguous arc of the circle, which is what makes the group
    colouring readable. The groups follow the order they are declared in
    config['parasite_groups'], so the circle and the legend read the same way.
    '''
    targets = {g: set(df['target']) for g, df in df_pred.groupby('taxid1_label')}
    labels = sorted(targets, key=lambda p: (group_order.get(groups.get(p), len(group_order)), p))

    nodes = [(i, g[0]+'. '+g.split(' ')[1], groups.get(g, UNKNOWN_GROUP)) for i, g in enumerate(labels)]
    links = [(i, j, len(targets[g1].intersection(targets[labels[j]])))
             for i, g1 in enumerate(labels) for j in range(i + 1, len(labels))]

    return (pd.DataFrame(links, columns=['source', 'target', 'value']),
            pd.DataFrame(nodes, columns=['index', 'name', 'group']))

@st.cache_resource(show_spinner=False)
def generate_circos_plot(data_dir, host_taxids, groups, palette, config, only_expressed=False):
    '''
    The arcs are coloured by taxonomic group (config['parasite_groups']) rather than
    by species: a per-species palette would have to repeat itself for the 35 parasites
    infecting human, and the group colour also shows whether the shared host proteins
    follow the parasite phylogeny. The chords themselves are neutral grey -- nearly
    every pair of parasites shares something, so colouring them only adds noise --
    and hovering (or clicking) a parasite draws the chords it belongs to in green.
    The hover styling has to set the alpha too: inheriting the 0.35 of the resting
    chords washes the highlight out to the same grey.

    With only_expressed the interactions are first restricted to the host proteins the
    parasite could meet in the tissues it infects.
    '''
    predictions = (get_tissue_expressed_predictions(data_dir, config, host_taxids) if only_expressed
                   else get_host_predictions(data_dir, host_taxids))
    links, nodes = get_common_interactors(predictions, groups,
                                          {g: i for i, g in enumerate(palette)})
    if links.empty or links['value'].max() < 1:
        return None, []

    palette = {**palette, **{g: UNKNOWN_COLOR for g in set(nodes['group']) if g not in palette}}
    # drop the pairs sharing nothing before scaling, so none of them is rounded up to a chord
    links = links[links['value'] >= 1].copy()
    links['value'] = spread_chord_strands(links['value'].to_numpy())

    chord = hv.Chord((links, hv.Dataset(nodes, 'index', ['name', 'group'])))
    chord.opts(
        opts.Chord(width=500, height=700, labels='name',
                   node_color=dim('group').str(), cmap=palette, node_line_color='white',
                   edge_line_color='#bdbdbd', edge_alpha=0.35,
                   edge_hover_line_color=HOVER_COLOR, edge_hover_line_alpha=1,
                   edge_selection_line_color=HOVER_COLOR, edge_selection_line_alpha=1,
                   edge_nonselection_line_alpha=0.1, inspection_policy='nodes'))

    shown = set(nodes['group'])

    figure = separate_arcs(inset_chord_ends(hv.render(chord)))

    return figure, [(g, c) for g, c in palette.items() if g in shown]

@st.cache_data(show_spinner=False)
def get_interactor_similarity(df_pred, groups):
    '''
    Jaccard similarity between the host proteins each pair of parasites is predicted to
    interact with, ordered by hierarchical clustering. Jaccard rather than the count of
    shared proteins the circos draws: the parasites have between 12 and 200 predicted
    interactors, so a raw count mostly reports how many predictions a parasite has.
    Clustering the parasites on that similarity puts the ones converging on the same host
    proteins next to each other, whether or not they are related.
    '''
    targets = {g: set(df['target']) for g, df in df_pred.groupby('taxid1_label')}
    targets = {g: t for g, t in targets.items() if t}
    labels = sorted(targets)
    if len(labels) < 3:
        return None

    similarity = np.array([[len(targets[a] & targets[b]) / len(targets[a] | targets[b])
                            for b in labels] for a in labels])
    order = leaves_list(linkage(squareform(1 - similarity, checks=False), method='average'))
    labels = [labels[i] for i in order]

    return (pd.DataFrame(similarity[np.ix_(order, order)], index=labels, columns=labels),
            [groups.get(g, UNKNOWN_GROUP) for g in labels])


@st.cache_data(show_spinner=False)
def generate_similarity_heatmap(similarity, clades, palette):
    '''
    The similarity matrix, with a strip of the taxonomic group of each parasite down the
    side so it can be read against the clustering. The diagonal is left blank: a parasite
    shares everything with itself and would take over the colour scale.
    '''
    shown = [g for g in palette if g in set(clades)]
    steps = [[i / len(shown), palette[g]] for i, g in enumerate(shown)]
    steps += [[(i + 1) / len(shown), palette[g]] for i, g in enumerate(shown)]
    names = [f'{g[0]}. {g.split(" ")[1]}' for g in similarity.index]

    figure = make_subplots(rows=1, cols=2, column_widths=[0.03, 0.97],
                           shared_yaxes=True, horizontal_spacing=0.01)
    figure.add_trace(go.Heatmap(z=[[shown.index(c)] for c in clades], y=names,
                                text=[[c] for c in clades], hovertemplate='%{text}<extra></extra>',
                                colorscale=sorted(steps), zmin=-0.5, zmax=len(shown) - 0.5,
                                showscale=False, xgap=0, ygap=1), row=1, col=1)

    values = similarity.to_numpy().copy()
    np.fill_diagonal(values, np.nan)
    figure.add_trace(go.Heatmap(z=values, x=names, y=names, colorscale='Blues', zmin=0, zmax=1,
                                hovertemplate='%{y} and %{x}<br>Jaccard %{z:.2f}<extra></extra>',
                                colorbar=dict(title='Shared<br>interactors<br>(Jaccard)',
                                              thickness=12, len=0.6)), row=1, col=2)

    # the strip is a heatmap and cannot carry a legend of its own, so the groups are named
    # by empty traces whose only purpose is their legend entry
    for group in shown:
        figure.add_trace(go.Scatter(x=[None], y=[None], mode='markers', name=group,
                                    marker=dict(size=10, symbol='square', color=palette[group]),
                                    hoverinfo='skip', showlegend=True), row=1, col=2)

    figure.update_xaxes(showticklabels=False, row=1, col=1)
    figure.update_xaxes(tickangle=-60, row=1, col=2)
    figure.update_yaxes(autorange='reversed')
    figure.update_layout(height=max(420, 22 * len(names) + 230), plot_bgcolor='white',
                         margin=dict(l=0, r=0, t=40, b=10),
                         legend=dict(orientation='h', yanchor='bottom', y=1.01, xanchor='left',
                                     x=0, itemclick=False, itemdoubleclick=False,
                                     font=dict(size=11)))

    return figure


def spread_chord_strands(values, target=480):
    '''
    Holoviews draws a chord as one line per unit of its value rather than as a filled
    ribbon, and spreads those lines evenly around the circle. Where the values are small --
    the hosts other than human, and any view of few parasites -- the lines end up degrees
    apart and a single chord reads as a handful of separate links. Scale the values up so
    that every host draws roughly the same number of lines, close enough together to close
    into a band. The proportions between the chords, and so the width of each, are kept.

    Holoviews scales down anything above its own max_chords of 500 itself, so the target
    stays below that and large values are left alone.
    '''
    total = values.sum()
    if total >= target:
        return values

    return np.maximum(1, np.round(values * target / total)).astype('int64')

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
        placed = [[] for _ in spans]
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
                    placed[node].append(angle)
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
        set_chord_width(renderer.edge_renderer, placed)

    return figure

def set_chord_width(edge_renderer, placed, radius=250):
    '''
    Widen each line of a chord until it meets its neighbour, so a chord reads as one band
    of a width rather than as a row of separate links. The spacing has to be measured from
    the drawn chords: holoviews scales the values down to its max_chords before deciding
    how many lines to draw, so it does not follow from the number of shared proteins.
    '''
    apart = [np.diff(np.sort(angles)) for angles in placed if len(angles) > 1]
    if not apart:
        return

    width = float(min(max(np.median(np.concatenate(apart)) * radius, 0.6), 3.0))
    edge_renderer.glyph.line_width = width
    for state in ('hover_glyph', 'selection_glyph'):
        glyph = getattr(edge_renderer, state, None)
        if glyph is not None:
            glyph.line_width = width * 1.6

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

def show_circos_plot(data_dir, host, host_taxids, config, only_expressed, caption, key):
    st.caption(caption)
    circos_plot, circos_legend = generate_circos_plot(
        data_dir, host_taxids,
        {p['label']: p.get('group', UNKNOWN_GROUP) for p in config['parasites'].values()},
        config.get('parasite_groups', {}), config, only_expressed)
    if circos_plot is None:
        where = ' in the tissues they infect' if only_expressed else ''
        st.text(f'No host proteins are shared by the parasites infecting {host}{where}')
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
        shared = (f'Each arc is a parasite infecting {selected_host}, coloured by its taxonomic '
                  'group; a chord joins two parasites that are predicted to interact with the '
                  'same host proteins, and is the thicker the more proteins they share. Hover a '
                  'parasite to follow its chords. ')
        views = {
            IN_TISSUE_VIEW: shared + 'Only the host proteins expressed in the tissues each '
                                     'parasite is known to infect are counted, so what is left '
                                     'are the interactions that can take place where the '
                                     'parasite is.',
            'All predictions': shared + 'Every predicted interaction is counted, whether or not '
                                        'the host protein is expressed where the parasite is.',
        }
        # one view at a time rather than two tabs: the bokeh component sizes itself from the
        # width of the document, which is 0 while it sits in a tab that is not open
        view = st.radio('Interactions to count', list(views), horizontal=True,
                        key='circos_view', label_visibility='collapsed')
        show_circos_plot(data_dir, selected_host, selected_taxids, config,
                         view == IN_TISSUE_VIEW, views[view], 'circos')

    # the heatmap answers the same question as the circos, on the interactions the radio
    # above selects, but reads where the circos cannot: it normalises for how many
    # predictions a parasite has, and it has no chords to overlap
    similarity = get_interactor_similarity(
        get_tissue_expressed_predictions(data_dir, config, selected_taxids)
        if view == IN_TISSUE_VIEW else get_host_predictions(data_dir, selected_taxids),
        {p['label']: p.get('group', UNKNOWN_GROUP) for p in config['parasites'].values()})

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
                   'many predictions does not look similar to everything. The parasites are '
                   'ordered by clustering them on that similarity and the strip on the left is '
                   'their taxonomic group, so a clade appearing as one block means its parasites '
                   'converge on the same host proteins.')
        st.plotly_chart(generate_similarity_heatmap(*similarity, config.get('parasite_groups', {})),
                        width='stretch')

st.markdown("---")


# Footer
with st.container():
    web_utils.footer()