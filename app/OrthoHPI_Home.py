import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils
import web_utils
import streamlit as st
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
def get_host_predictions(data_dir, host):
    predictions = load_predictions(data_dir)

    return predictions[predictions['taxid2_label'] == host]


@st.cache_data(show_spinner=False)
def count_interactions_per_tissue_cell_type(data_dir, config, host):
    '''
    Counts the predicted interactions of a host per parasite, tissue and cell type,
    keeping only the tissues each parasite is known to infect (config['parasites']).
    Aggregating here keeps the merged predictions/tissues table out of the app's
    memory and makes sure each icicle sector is counted once.
    '''
    predictions = get_host_predictions(data_dir, host)
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
def get_common_interactors(df_pred):
    targets = {g: set(df['target']) for g, df in df_pred.groupby('taxid1_label')}
    labels = list(targets)

    nodes = [(i, g[0]+'. '+g.split(' ')[1]) for i, g in enumerate(labels)]
    links = [(i, j, len(targets[g1].intersection(targets[labels[j]])))
             for i, g1 in enumerate(labels) for j in range(i + 1, len(labels))]

    return (pd.DataFrame(links, columns=['source', 'target', 'value']),
            pd.DataFrame(nodes, columns=['index', 'name']))

@st.cache_resource(show_spinner=False)
def generate_circos_plot(data_dir, host):
    links, nodes = get_common_interactors(get_host_predictions(data_dir, host))
    if links.empty or links['value'].max() < 1:
        return None

    chord = hv.Chord((links, hv.Dataset(nodes, 'index'))).select(value=(1, None))
    chord.opts(
        opts.Chord(width=500, height=700, cmap='Category20', edge_cmap='Category20', edge_color=dim('source').str(),
               labels='name', node_color=dim('index').str()))

    return hv.render(chord)

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
    stats_figures.append((fig, "Number of Interactions"))

    fig = generate_boxplot_score_stats(df)
    stats_figures.append((fig, "Boxplot of Confidence scores"))

    return stats_figures




st.markdown("<h1 style='text-align: center; color: #023858;'>OrthoHPI 2.0</h1>", unsafe_allow_html=True)
st.markdown("<h3 style='text-align: center; color: #2b8cbe;'>Orthology Prediction of Host-Parasite PPI</h3>", unsafe_allow_html=True)

st.text(" ")
st.text(" ")
st.markdown("---")

# Define selection options
host_list = ['<select>'] + load_predictions(data_dir)['taxid2_label'].sort_values().unique().tolist()

col1, col2, col3 = st.columns(3)

with col1:
    st.write('')

with col2:
    selected_host = st.selectbox('Select a host to visualize the predicted interactions', host_list, key="home_host")
    if selected_host == "<select>":
        st.text('Choose 1 host to explore the predicted host-parasite interactions')

with col3:
    st.write('')


if selected_host != "<select>":
    chart1, chart2 = st.columns(2)

    with chart1:
        st.subheader("Circos Plot of Common Host Interactors")
        circos_plot = generate_circos_plot(data_dir, selected_host)
        if circos_plot is not None:
            streamlit_bokeh(circos_plot, use_container_width=True)
        else:
            st.text(f'No host proteins are shared by the parasites infecting {selected_host}')

    stats_figs = generate_stats_plots(get_host_predictions(data_dir, selected_host))
    stats_cols = st.columns(len(stats_figs))
    i = 0
    for stats_fig, title in stats_figs:
        with stats_cols[i]:
            st.subheader(title)
            st.plotly_chart(stats_fig, width='stretch')
        i += 1

    fig = generate_tissue_cell_type_box(count_interactions_per_tissue_cell_type(data_dir, config, selected_host))
    with chart2:
        st.subheader("Summary of Interactions per Tissue and Cell type")
        st.plotly_chart(fig, width='stretch')

st.markdown("---")


# Footer
with st.container():
    web_utils.footer()