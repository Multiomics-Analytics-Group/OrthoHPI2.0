import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils
import web_utils
import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
import holoviews as hv
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

    chord = hv.Chord((links, hv.Dataset(nodes, 'index', ['name', 'group']))).select(value=(1, None))
    chord.opts(
        opts.Chord(width=500, height=700, labels='name',
                   node_color=dim('group').str(), cmap=palette, node_line_color='white',
                   edge_line_color='#bdbdbd', edge_line_width=0.6, edge_alpha=0.35,
                   edge_hover_line_color=HOVER_COLOR, edge_hover_line_alpha=1,
                   edge_hover_line_width=1.4,
                   edge_selection_line_color=HOVER_COLOR, edge_selection_line_alpha=1,
                   edge_selection_line_width=1.4, edge_nonselection_line_alpha=0.1,
                   inspection_policy='nodes'))

    shown = set(nodes['group'])

    return separate_arcs(hv.render(chord)), [(g, c) for g, c in palette.items() if g in shown]

def separate_arcs(figure, gap=0.008):
    '''
    Holoviews draws each node arc from where the previous one ends, with no gap. Now that
    neighbouring parasites of the same taxonomic group share a colour, they merge into one
    band, and a chord landing on either of two neighbours looks like two chords to the same
    parasite. Trim every arc by a fixed angle, clamped so short arcs are shortened rather
    than erased, so each parasite reads as its own segment.
    '''
    trimmed = set()
    for renderer in figure.renderers:
        data = getattr(renderer, 'data_source', None) and renderer.data_source.data
        if not data or 'arc_xs' not in data or id(data) in trimmed:
            continue
        trimmed.add(id(data))

        arc_xs, arc_ys = [], []
        for xs, ys in zip(data['arc_xs'], data['arc_ys']):
            angles = np.unwrap(np.arctan2(ys, xs))
            extent = angles[-1] - angles[0]
            pad = np.sign(extent) * min(gap, abs(extent) * 0.3)
            angles = np.linspace(angles[0] + pad, angles[-1] - pad, len(angles))
            arc_xs.append(np.cos(angles))
            arc_ys.append(np.sin(angles))
        data['arc_xs'], data['arc_ys'] = arc_xs, arc_ys

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

st.markdown("---")


# Footer
with st.container():
    web_utils.footer()