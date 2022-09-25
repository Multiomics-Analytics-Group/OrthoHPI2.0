import os
import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import networkx as nx
from pyvis.network import Network

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
        value = v*2/0.9
        if value < 0.2:
            value = 0.2
        widths[(n1, n2)] = value
    nx.set_edge_attributes(G, widths, 'value')
    

    return G

# Read dataset
predictions = pd.read_csv('data/predictions.tsv', sep='\t', header=0)
tissues = pd.read_csv('data/tissues_cell_types.tsv', sep='\t', header=0)


st.title('Human-parsite Protein-Protein Interaction Graph')

# Define selection options and sort alphabetically
parasite_list = ['<select>'] + predictions['taxid1_label'].unique().tolist()

# Implement multiselect dropdown menu for option selection
selected_parasite = st.selectbox('Select a parasite to visualize the predicted PPI', parasite_list)

# Set info message on initial site load
if selected_parasite == "<select>":
   st.text('Please choose at least 1 parasite to get started')
else:
    score = st.slider('Confidence score', 0.4, 0.9, 0.7)
    
    df_select = predictions.loc[predictions['taxid1_label'] == selected_parasite]  
    # Create networkx graph object from pandas dataframe
    G = generate_graph(df_select, score)

    st.text(f"Nodes: {len(G.nodes())}  Edges: {len(G.edges())}")

    # Initiate PyVis network object
    net = Network(height='1000px', width='1000px', bgcolor='#222222', font_color='white')

    # Take Networkx graph and translate it to a PyVis graph format
    net.from_nx(G)

    # Generate network with specific layout settings
    net.repulsion(node_distance=420, central_gravity=0.33,
                       spring_length=110, spring_strength=0.10,
                       damping=0.95)
    
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
    components.html(HtmlFile.read(), height=1000, width=1000)

