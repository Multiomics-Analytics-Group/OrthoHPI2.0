import os
import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import networkx as nx
from pyvis.network import Network
import plotly.express as px
import holoviews as hv
from holoviews import opts, dim
hv.extension('bokeh')


st.set_page_config(layout="wide")
st.markdown("<h1 style='text-align: center; color: #023858;'>OrthoHPI 2.0</h1>", unsafe_allow_html=True)
st.markdown("<h3 style='text-align: center; color: #2b8cbe;'>Homology Prediction of Host-Parasite PPI</h3>", unsafe_allow_html=True)

@st.cache(suppress_st_warning=True, allow_output_mutation=True)
def generate_graph(df, score):
    G = nx.from_pandas_edgelist(df, 'source_name', 'target_name', 'weight')
    colors = dict(df[['source_name', 'source_color']].drop_duplicates().values)
    colors.update(dict(df[['target_name', 'target_color']].drop_duplicates().values))
    nx.set_node_attributes(G, colors, 'color')
    shapes = dict(df[['source_name', 'source_shape']].drop_duplicates().values)
    shapes.update(dict(df[['target_name', 'target_shape']].drop_duplicates().values))
    nx.set_node_attributes(G, shapes, 'shape')
    centrality = nx.betweenness_centrality(G, weight='weight')
    max_centrality = max(list(centrality.values()))
    sizes = {}
    for k,v in centrality.items():
        value = v*60/max_centrality
        if value < 20:
            value = 20
        sizes[k] =  value
    nx.set_node_attributes(G, sizes, 'size')

    rm_edges = [(n1, n2) for n1,n2,w in G.edges.data('weight') if w < score]
    # remove filtered edges from graph G
    G.remove_edges_from(rm_edges)
    G.remove_nodes_from(list(nx.isolates(G)))

    widths = {}
    for n1,n2,w in df[['source_name', 'target_name', 'weight']].values:
        value = v*0.5/0.9
        if value < 0.05:
            value = 0.05
        widths[(n1, n2)] = value
    nx.set_edge_attributes(G, widths, 'value')
    nx.set_edge_attributes(G, '#999999', 'color')
    

    return G

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


def generate_tissue_filters(df):
    options = df['Tissue'].unique().tolist()
    
    return options

def generate_cell_type_filters(df):
    options = df['Cell type'].dropna().unique().tolist()
    
    return options

# Read dataset
predictions = pd.read_csv('data/predictions.tsv', sep='\t', header=0)
tissues = pd.read_csv('data/tissues_cell_types.tsv', sep='\t', header=0)
pred_tissues = pd.merge(predictions.rename({'target':'Gene'}, axis=1), tissues, on='Gene', how='left')

df_select = None
net = None

# Define selection options
parasite_list = ['<select>'] + predictions['taxid1_label'].unique().tolist()

chart1, chart2 = st.columns(2)


with chart1:
    st.subheader("Circos Plot of Common Host Interactors")
    circos_plot = generate_circos_plot(predictions)
    st.bokeh_chart(hv.render(circos_plot, backend='bokeh'), use_container_width=True)

stats_figs = generate_stats_plots(predictions)
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



st.markdown("<h3 style='text-align: center; color: black;'>Graph of predicted Host-Parasite PPIs</h3>", unsafe_allow_html=True)


col1, col2, col3 = st.columns(3)

with col1:
    st.write('')


with col2:
    
    # Implement multiselect dropdown menu for option selection
    selected_parasite = st.selectbox('Select a parasite to visualize the predicted PPI', parasite_list)

    # Set info message on initial site load
    if selected_parasite == "<select>":
        st.text('Choose 1 parasitic specie to visualize the predicted PPI network')
    else:        
        df_select = pred_tissues.loc[pred_tissues['taxid1_label'] == selected_parasite]
        score = st.slider('Confidence score', 0.4, 0.9, 0.7)
        
        tissues_options = generate_tissue_filters(pred_tissues)
        if len(tissues_options) > 0:
            selected_tissues = st.multiselect('Select tissues to filter the predicted PPI', tissues_options)
            if len(selected_tissues) > 0:
                df_select = df_select[df_select['Tissue'].isin(selected_tissues)]
                cell_type_options = generate_cell_type_filters(df_select)
                if len(cell_type_options) > 0:
                    selected_cell_types = st.multiselect('Select cell type to filter the predicted PPI', cell_type_options)
                    if len(selected_cell_types) > 0 :
                        df_select = df_select[df_select['Cell type'].isin(selected_cell_types)]


        # Create networkx graph object from pandas dataframe
        G = generate_graph(df_select, score)

        st.text(f"Nodes: {len(G.nodes())}  Edges: {len(G.edges())}")

        # Initiate PyVis network object
        net = Network(height='1000px', width="100%", bgcolor='white', font_color='#555555')
        # Take Networkx graph and translate it to a PyVis graph format
        net.from_nx(G)

        # Generate network with specific layout settings
        net.repulsion(node_distance=420, central_gravity=0.33,
                        spring_length=110, spring_strength=0.10,
                        damping=0.95)
        
        #net.show_buttons(filter_=['nodes'])
        
        
with col3:
    st.write('')


with st.container():
    if net is not None:
        # Save and read graph as HTML file (on Streamlit Sharing)
        try:
            path = 'data/tmp'
            net.save_graph(f'{path}/{selected_parasite}.html')
            HtmlFile = open(f'{path}/{selected_parasite}.html','r',encoding='utf-8')
            # Save and read graph as HTML file (locally)
        except:
            path = 'data/html_files'
            net.save_graph(f'{path}/{selected_parasite}.html')
            HtmlFile = open(f'{path}/{selected_parasite}.html','r',encoding='utf-8')

        # Load HTML into HTML component for display on Streamlit
        components.html(HtmlFile.read(), height=1200)

with st.container():
    if df_select is not None:
        st.header("Table of Host-Parasite PPIs")
        st.dataframe(df_select)
