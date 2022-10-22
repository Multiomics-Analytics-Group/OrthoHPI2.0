import utils
import web_utils
import streamlit as st
import pandas as pd
import plotly.express as px
import holoviews as hv
from css import style
from holoviews import opts, dim
hv.extension('bokeh')

st.set_page_config(layout="wide", page_title="OrthoHPI Home", menu_items={})

style.load_css()


# Read dataset
config = utils.read_config('config.yml')
predictions = utils.read_parquet_file(input_file='data/predictions.parquet')
predictions['weight'] = predictions['weight'].astype(float)
tissues = utils.read_parquet_file(input_file='data/tissues_cell_types.parquet')
pred_tissues = pd.merge(predictions, tissues.rename({'Gene': 'target'}, axis=1), on='target', how='left')
tissues = None
ontology = utils.read_parquet_file(input_file='data/go_ontology.parquet')

#Initialize variables
df_select = None
net = None
selected_rows = []
selected_terms = []
enrichment_table = None
enrichment = None


@st.cache(suppress_st_warning=True, allow_output_mutation=True)
def generate_tissue_cell_type_box(df):
    aux = df.copy()
    aux['Cell type'] = aux['Cell type'].fillna("Not available")
    counts_tissues = aux.groupby(['taxid1', 'Tissue']).count()['taxid2'].reset_index()
    counts_tissues = counts_tissues.rename({'taxid2':'edges_tissue'}, axis=1)
    counts_cells = aux.groupby(['taxid1', 'Tissue', 'Cell type']).count()['taxid2'].reset_index()
    counts_cells = counts_cells.rename({'taxid2':'edges_cell_type'}, axis=1)
    aux = pd.merge(aux, counts_tissues, on=['taxid1', 'Tissue'], how='left')
    aux = pd.merge(aux, counts_cells, on=['taxid1', 'Tissue', 'Cell type'], how='left')
    fig = px.icicle(aux, path=[px.Constant("Parasites"), 'taxid1_label', 'Tissue', 'Cell type'], values='edges_cell_type',
                  color='edges_cell_type', hover_data=['edges_tissue', 'edges_cell_type', 'taxid1', 'taxid1_label', 'pTPM'],
                  color_continuous_scale='Burgyl', height=900, width=1200)

    return fig

def generate_circos_plot(df_pred):
    nodes = set()
    links = []
    seen = set()
    i = 0
    for g1, df1 in df_pred.groupby('taxid1_label'):
        j = i + 1
        for g2, df2 in df_pred.groupby('taxid1_label'):
            if g1 != g2 and (g1, g2) not in seen:
                nodes.update([(i, g1[0]+'. '+g1.split(' ')[1]), (j, g2[0]+'. '+g2.split(' ')[1])])
                links.append((i, j, len(set(df1['target'].tolist()).intersection(df2['target'].tolist()))))
                seen.update([(g1, g2), (g2, g1)])
                
                j += 1
        i += 1

    links = pd.DataFrame(links, columns=['source', 'target', 'value'])
    nodes = hv.Dataset(pd.DataFrame(list(nodes), columns = ['index', 'name']), 'index')

    chord = hv.Chord((links, nodes)).select(value=(1, None))
    chord.opts(
        opts.Chord(width=500, height=700, cmap='Category20', edge_cmap='Category20', edge_color=dim('source').str(), 
               labels='name', node_color=dim('index').str()))

    return chord

def generate_boxplot_score_stats(df):
    fig = px.box(df.sort_values("taxid1"), x="taxid1_label", y="weight", color='taxid1', labels={"weight":"score", "taxid1_label": "parasites"})
    fig.update_traces(showlegend=False)

    return fig

def generate_barplot_stats(df):
    fig = px.bar(df.groupby(["taxid1_label"]).count().reset_index().sort_values("taxid2"), x="taxid1_label", y="weight", color='taxid1_label',  labels={"weight":"count", "taxid1_label": "parasites"})
    fig.update_traces(showlegend=False)
    return fig

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

chart1, chart2 = st.columns(2)


with chart1:
    st.subheader("Circos Plot of Common Host Interactors")
    circos_plot = generate_circos_plot(predictions)
    st.bokeh_chart(hv.render(circos_plot, backend='bokeh'), use_container_width=True)

stats_figs = generate_stats_plots(predictions)
predictions = None
stats_cols = st.columns(len(stats_figs))
i = 0
for stats_fig, title in stats_figs:
    with stats_cols[i]:
        st.subheader(title)
        st.plotly_chart(stats_fig, use_container_width=True)
    i += 1

fig = generate_tissue_cell_type_box(pred_tissues)
with chart2:
    st.subheader("Summary of Interactions per Tissue and Cell type")
    st.plotly_chart(fig, use_container_width=True)
    pred_tissues = None

st.markdown("---")


# Footer
with st.container():
    web_utils.footer()