import utils
import web_utils
import streamlit as st
from streamlit_extras.switch_page_button import switch_page
import streamlit.components.v1 as components
from st_aggrid import GridOptionsBuilder, AgGrid
import pandas as pd
import networkx as nx
from css import style
from pyvis.network import Network
import plotly.express as px

style.load_css()
page = web_utils.show_pages_menu(index=1)
if page == "Home":
    switch_page("orthohpi home")
elif page == "Predicted PPI structures":
    switch_page('interaction structures')
elif page == "About":
    switch_page('about')


#Initialize variables
df_select = None
net = None
selected_rows = []
selected_terms = []
enrichment_table = None
enrichment = None
path = 'data/tmp'

# Read dataset
config = utils.read_config('config.yml')
predictions = utils.read_parquet_file(input_file='data/predictions.parquet')
predictions['weight'] = predictions['weight'].astype(float)
tissues = utils.read_parquet_file(input_file='data/tissues_cell_types.parquet')
pred_tissues = pd.merge(predictions, tissues.rename({'Gene': 'target'}, axis=1), on='target', how='left')
tissues = None
predictions = None
ontology = utils.read_parquet_file(input_file='data/go_ontology.parquet')


def generate_tissue_filters(df):
    options = df['Tissue'].unique().tolist()
    
    return options

def generate_cell_type_filters(df):
    options = df['Cell type'].dropna().unique().tolist()
    
    return options

@st.cache_data
def get_enrichment(pred_df):
    species = pred_df['taxid1'].unique().tolist() + pred_df['taxid2'].unique().tolist()
    species = [int(s) for s in species]
    go_df = utils.read_parquet_file(input_file='data/gos.parquet')
    go_df = go_df[go_df['taxid'].isin(species)]
    enrichment = utils.calculate_enrichment(pred_df, go_df)

    return enrichment


def get_enrichment_summary(enrichment_df, ontology_df):
    df = ontology_df[(ontology_df['parent'].isin(enrichment_df['go_term'])) & (ontology_df['child'].isin(enrichment_df['go_term']))]
    df = pd.merge(df.rename({'child':'go_term'}, axis=1), enrichment_df[['go_term', 'odds_ratio', 'fdr_bh']], on='go_term')
    fig = px.treemap(df, path=['parent', 'go_term'], values='odds_ratio', height=900, hover_data=['fdr_bh', 'odds_ratio'])

    return fig

def generate_graph(df, score):
    G = nx.from_pandas_edgelist(df, 'source', 'target', 'weight')
    colors = dict(df[['source', 'source_color']].drop_duplicates().values)
    colors.update(dict(df[['target', 'target_color']].drop_duplicates().values))
    nx.set_node_attributes(G, colors, 'color')
    labels = dict(df[['source', 'source_name']].drop_duplicates().values)
    labels.update(dict(df[['target', 'target_name']].drop_duplicates().values))
    nx.set_node_attributes(G, labels, 'label')
    shapes = dict(df[['source', 'source_shape']].drop_duplicates().values)
    shapes.update(dict(df[['target', 'target_shape']].drop_duplicates().values))
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
    for n1,n2,w in df[['source', 'target', 'weight']].values:
        value = v*0.5/0.9
        if value < 0.05:
            value = 0.05
        widths[(n1, n2)] = value
    nx.set_edge_attributes(G, widths, 'value')
    nx.set_edge_attributes(G, '#999999', 'color')
    

    return G


st.markdown("<h1 style='text-align: center; color: #023858;'>OrthoHPI 2.0</h1>", unsafe_allow_html=True)
st.markdown("<h3 style='text-align: center; color: #2b8cbe;'>Orthology Prediction of Host-Parasite PPI</h3>", unsafe_allow_html=True)


# Define selection options
parasite_list = ['<select>'] + pred_tissues['taxid1_label'].sort_values().unique().tolist()

st.markdown("<h3 style='text-align: center; color: black;'>Graph of predicted Host-Parasite PPIs</h3>", unsafe_allow_html=True)


col1, col2, col3 = st.columns(3)

with col1:
    st.write('')


with col2:
    
    # Implement multiselect dropdown menu for option selection
    selected_parasite = st.selectbox('Select a parasite to visualize the predicted PPI', parasite_list, key="net_par")

    # Set info message on initial site load
    if selected_parasite == "<select>":
        st.text('Choose 1 parasite to visualize the predicted PPI network')
    else:        
        df_select = pred_tissues.loc[pred_tissues['taxid1_label'] == selected_parasite]
        pred_tissues = None
        df_select = web_utils.filter_tissues(config, df_select)
        score = st.slider('Confidence score', 0.4, 0.9, 0.7)

        tissues_options = generate_tissue_filters(df_select)
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
        # Save other formats
        utils.export_graph(G, filename=f'{selected_parasite}.graphml',
                        format='graphml', output_dir=f'{path}')
        utils.export_graph(G, filename=f'{selected_parasite}.json',
                        format='cytoscape', output_dir=f'{path}')
        G = None

        # Generate network with specific layout settings
        net.repulsion(node_distance=420, central_gravity=0.33,
                        spring_length=110, spring_strength=0.10,
                        damping=0.95)
        
        #net.show_buttons(filter_=['nodes'])
        
        
with col3:
    st.write('')


with st.container():
    if net is not None:
        html_data = ""
        # Save and read graph as HTML file (on Streamlit Sharing)
        net.save_graph(f'{path}/{selected_parasite}.html')
        with open(f'{path}/{selected_parasite}.html','r',encoding='utf-8') as HtmlFile:
            html_data = HtmlFile.read()
        # Load HTML into HTML component for display on Streamlit
        components.html(html_data, height=1050)
        net = None
        with st.container():
            c1, c2, c3 = st.columns(3)

            with c1:
                st.download_button(
                    label="Download Network as Html",
                    data=html_data,
                    file_name=f'{selected_parasite}_network.html',
                    mime='text/html',
                )
            with c2:
                st.download_button(
                    label="Download Network as GraphML",
                    data=open(f'{path}/{selected_parasite}.graphml','r',encoding='utf-8'),
                    file_name=f'{selected_parasite}_network.graphml',
                    mime='text/plain',
                )
            with c3:
                st.download_button(
                    label="Download Network as Cytoscape",
                    data=open(f'{path}/{selected_parasite}.json','r',encoding='utf-8'),
                    file_name=f'{selected_parasite}_network.json',
                    mime='text/plain',
                )


with st.container():
    if df_select is not None:
        st.header("Table of Host-Parasite PPIs")
        table = df_select[df_select['weight'] >= score]
        table = table.sort_values(by='weight', ascending=False)
        gb = GridOptionsBuilder.from_dataframe(table)
        gb.configure_pagination(paginationAutoPageSize=True) #Add pagination
        gb.configure_side_bar() #Add a sidebar
        gridOptions = gb.build()
        grid_response = AgGrid(
                            table,
                            gridOptions=gridOptions,
                            data_return_mode='AS_INPUT', 
                            fit_columns_on_grid_load=False,
                            enable_enterprise_modules=True,
                            height=350, 
                            reload_data=False
                        )
        st.download_button(
            label="Download Network Table",
            data=utils.convert_df(table),
            file_name=f'{selected_parasite}_network_table.tsv',
            mime='text/csv',
        )

with st.container():
    if df_select is not None:
        st.header("Network Functional Enrichment -- GO Biological Processes")
        enrichment = get_enrichment(df_select[df_select['weight'] >= score])
        if not enrichment.empty:
            fdr = st.radio("FDR BH correction",(0.01, 0.05, 0.1), horizontal=True)
            st.text(f"Terms enriched: {len(enrichment[enrichment['fdr_bh'] <= fdr]['go_term'].values.tolist())}")
            st.text("Select GO terms to get more details")
            enrichment_table = enrichment[enrichment['fdr_bh'] <= fdr][['go_term', 'p_value', 'odds_ratio', 'fdr_bh']]
            gb = GridOptionsBuilder.from_dataframe(enrichment_table)
            gb.configure_pagination(paginationAutoPageSize=True) #Add pagination
            gb.configure_side_bar() #Add a sidebar
            gb.configure_selection('multiple', use_checkbox=True, groupSelectsChildren="Group checkbox select children") #Enable multi-row selection
            gridOptions = gb.build()
            grid_response = AgGrid(
                                enrichment_table,
                                gridOptions=gridOptions,
                                data_return_mode='AS_INPUT', 
                                update_mode='MODEL_CHANGED', 
                                fit_columns_on_grid_load=False,
                                enable_enterprise_modules=True,
                                height=350, 
                                reload_data=False
                            )
            selected_rows = grid_response['selected_rows']
            st.download_button(
                label="Download Enrichment Table",
                data=utils.convert_df(enrichment_table),
                file_name=f'{selected_parasite}_network_enrichment_table.tsv',
                mime='text/csv',
            )
        else:
            st.subheader("No GO terms where found enriched")

go1, go2 = st.columns(2)
with st.container():
    if enrichment_table is not None:
        enrichment_viz = enrichment_table
        if selected_rows is not None and len(selected_rows) > 0:
            selected_terms = [i['go_term'] for i in selected_rows]
            enrichment_viz = enrichment_viz[enrichment_viz['go_term'].isin(selected_terms)]

        with go1:
            fig = px.scatter(enrichment_viz, x='fdr_bh', y='odds_ratio', 
                size='odds_ratio', color='go_term', height=450, 
                labels = {'fdr_bh':'FDR BH', 'odds_ratio': 'Odds ratio'})
            fig.update_traces(showlegend=False)
            st.subheader("Enriched Biological Processes -- Odds ratio vs FDR")
            st.plotly_chart(fig, height=400, use_container_width=True)
        with go2:
            if len(selected_terms) > 0:
                if enrichment is not None:
                    highlighted_nodes = enrichment[enrichment['go_term'].isin(selected_terms)]['nodes'].values
                    highlighted_nodes = utils.merge_list_of_lists([i.split(',') for i in highlighted_nodes])
                    highlight_color = {i: '#e7298a' for i in highlighted_nodes}
                    G = generate_graph(df_select, score)
                    nx.set_node_attributes(G, "#ddd", 'color')
                    nx.set_node_attributes(G, highlight_color, 'color')
                    # Initiate PyVis network object
                    net = Network(height="450px", width="100%", bgcolor='white', font_color='#555555')
                    # Take Networkx graph and translate it to a PyVis graph format
                    net.from_nx(G)
                    G = None
                    net.save_graph(f'{path}/{selected_parasite}2.html')
                    net = None
                    st.subheader("Highlighted Nodes for Selected Biological Processes")
                    with open(f'{path}/{selected_parasite}2.html','r',encoding='utf-8') as HtmlFile:
                        html_data = HtmlFile.read()
                    components.html(html_data, height=500)
                    st.download_button(
                        label="Download Network as Html",
                        data=html_data,
                        file_name=f'{selected_parasite}_enrichment_network.html',
                        mime='text/html',
                    )
        
        fig = get_enrichment_summary(enrichment_table, ontology)
        st.subheader("Visual Summary of Enriched Hierarchy of Biological Processes")
        st.plotly_chart(fig, use_container_width=True)

st.markdown("---")
st.markdown("---")

# Footer
with st.container():
    web_utils.footer()