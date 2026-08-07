import sys, os
import textwrap
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import utils
import web_utils
import streamlit as st
from st_aggrid import GridOptionsBuilder, AgGrid
import pandas as pd
import networkx as nx
from css import style
from pyvis.network import Network
import plotly.express as px

style.load_css()
page = web_utils.show_pages_menu(index=1)
if page == "Home":
    st.switch_page("OrthoHPI_Home.py")
elif page == "Predicted PPI structures":
    st.switch_page('pages/2_Interaction_structures.py')
elif page == "About":
    st.switch_page('pages/3_About.py')

#Initialize variables
df_select = None
net = None
selected_rows = []
selected_terms = []
enrichment_table = None
enrichment = None
path = 'data/tmp'
# characters per line of a node label, so a long protein name wraps instead of
# stretching its node across the network
LABEL_WRAP_WIDTH = 22
LABEL_FONT_SIZE = 22  # vis.js defaults to 14, too small to read the protein names
LABEL_FONT_COLOR = '#555555'

# Read dataset
config = utils.read_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()


@st.cache_data(show_spinner=False)
def load_predictions(data_dir):
    predictions = utils.read_parquet_file(input_file=f'{data_dir}/predictions.parquet')
    predictions['weight'] = predictions['weight'].astype(float)

    return predictions


@st.cache_data(show_spinner=False, max_entries=5)
def get_parasite_tissues(data_dir, parasite, host_taxids):
    '''
    Annotates the predictions of one parasite against the selected host with the
    tissues and cell types its host targets are expressed in. Selecting the parasite
    before merging keeps the full predictions x tissues table (~630 MB) out of the app.
    '''
    predictions = load_predictions(data_dir)
    tissues = utils.read_parquet_file(input_file=f'{data_dir}/tissues_cell_types.parquet')
    predictions = predictions[(predictions['taxid1_label'] == parasite)
                              & (predictions['taxid2'].isin(host_taxids))]

    return pd.merge(predictions, tissues.rename({'Gene': 'target'}, axis=1), on='target', how='left')


@st.cache_data(show_spinner=False)
def get_parasite_list(data_dir, host_taxids):
    '''Parasites with predictions against the selected host, so the parasite
    selector never offers one that yields an empty network.'''
    predictions = load_predictions(data_dir)
    predictions = predictions[predictions['taxid2'].isin(host_taxids)]

    return predictions['taxid1_label'].sort_values().unique().tolist()


@st.cache_data(show_spinner=False)
def load_ontology(data_dir):
    return utils.read_parquet_file(input_file=f'{data_dir}/go_ontology.parquet')


def set_label_font(net):
    '''
    Enlarges the node labels of a pyvis network. The size has to be set on the
    network rather than on each node: pyvis overwrites a node's font with the
    font_color the network was built with.

    :param net: pyvis Network whose labels to enlarge
    '''
    net.options.nodes = {'font': {'size': LABEL_FONT_SIZE, 'color': LABEL_FONT_COLOR}}


def generate_node_labels(df, annotations):
    '''
    Draws the descriptive protein name on the node instead of the short name STRING
    prefers. Proteins STRING has no name for are left with their short name -- an
    "Uncharacterized protein" label identifies nothing, and most parasite proteins
    would end up sharing it. The name is wrapped so long ones do not stretch the node
    across the network.

    :param df: predictions dataframe of the selected parasite
    :param dict annotations: STRING id --> descriptive protein name
    :return: {STRING id: label}
    '''
    labels = {}
    for prefix in ['source', 'target']:
        for protein, name in df[[prefix, f'{prefix}_name']].drop_duplicates(subset=prefix).values:
            description = annotations.get(protein, '')
            if not description or description.lower().startswith('uncharacterized'):
                labels[protein] = str(name)
            else:
                labels[protein] = '\n'.join(textwrap.wrap(description, width=LABEL_WRAP_WIDTH))

    return labels


def generate_node_titles(df, annotations):
    '''
    Builds the hover text of each node: the short name, the descriptive protein name
    and the identifiers needed to look the protein up elsewhere. The label only shows
    one of the two names, so the tooltip keeps both.

    :param df: predictions dataframe of the selected parasite
    :param dict annotations: STRING id --> descriptive protein name
    :return: {STRING id: hover text}
    '''
    titles = {}
    for prefix, taxid_col in [('source', 'taxid1_label'), ('target', 'taxid2_label')]:
        cols = [prefix, f'{prefix}_name', f'{prefix}_uniprot', taxid_col]
        for protein, name, uniprot, species in df[cols].drop_duplicates(subset=prefix).values:
            lines = [str(name)]
            description = annotations.get(protein)
            if description and description != name:
                lines.append(description)
            lines.append(str(species))
            lines.append(f'STRING: {protein}')
            if pd.notna(uniprot):
                lines.append(f'UniProt: {uniprot}')
            titles[protein] = '\n'.join(lines)

    return titles


def generate_tissue_filters(df):
    options = df['Tissue'].unique().tolist()
    
    return options

def generate_cell_type_filters(df):
    options = df['Cell type'].dropna().unique().tolist()
    
    return options

@st.cache_data(max_entries=3, ttl=1800)
def get_enrichment(pred_df, data_dir):
    species = pred_df['taxid1'].unique().tolist() + pred_df['taxid2'].unique().tolist()
    species = [int(s) for s in species]
    # the filter is pushed down to the reader (fastparquet prunes row groups only,
    # so the exact selection is still applied afterwards)
    go_df = utils.read_parquet_file(input_file=f'{data_dir}/gos.parquet', filters=[('taxid', 'in', species)])
    go_df = go_df[go_df['taxid'].isin(species)]
    enrichment = utils.calculate_enrichment(pred_df, go_df)

    return enrichment


def get_enrichment_summary(enrichment_df, ontology_df):
    df = ontology_df[(ontology_df['parent'].isin(enrichment_df['go_term'])) & (ontology_df['child'].isin(enrichment_df['go_term']))]
    df = pd.merge(df.rename({'child':'go_term'}, axis=1), enrichment_df[['go_term', 'odds_ratio', 'fdr_bh']], on='go_term')
    fig = px.treemap(df, path=['parent', 'go_term'], values='odds_ratio', height=900, hover_data=['fdr_bh', 'odds_ratio'])

    return fig

def generate_graph(df, score, annotations=None):
    if df.empty:
        return nx.Graph()

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
    if annotations is not None:
        nx.set_node_attributes(G, generate_node_labels(df, annotations), 'label')
        nx.set_node_attributes(G, generate_node_titles(df, annotations), 'title')
    centrality = nx.betweenness_centrality(G, weight='weight')
    max_centrality = max(centrality.values(), default=0)
    sizes = {}
    for k,v in centrality.items():
        value = v*60/max_centrality if max_centrality > 0 else 20
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
        value = w*0.5/0.9
        if value < 0.05:
            value = 0.05
        widths[(n1, n2)] = value
    nx.set_edge_attributes(G, widths, 'value')
    nx.set_edge_attributes(G, '#999999', 'color')
    

    return G


st.markdown("<h1 style='text-align: center; color: #023858;'>OrthoHPI 2.0</h1>", unsafe_allow_html=True)
st.markdown("<h3 style='text-align: center; color: #2b8cbe;'>Orthology Prediction of Host-Parasite PPI</h3>", unsafe_allow_html=True)


st.markdown("<h3 style='text-align: center; color: black;'>Graph of predicted Host-Parasite PPIs</h3>", unsafe_allow_html=True)


col1, col2, col3 = st.columns(3)

with col1:
    st.write('')


with col2:

    # the host carries over from whichever page it was last chosen on
    selected_host, selected_taxids = web_utils.host_selector(
        config, load_predictions(data_dir), 'Select a host')

    if selected_host == web_utils.NO_HOST:
        st.text('Choose 1 host to explore the predicted host-parasite interactions')
        selected_parasite = "<select>"
    else:
        # only parasites that infect the selected host
        parasite_list = ['<select>'] + get_parasite_list(data_dir, selected_taxids)
        # switching host can leave a parasite selected that the new host does not have
        if st.session_state.get('net_par') not in parasite_list:
            st.session_state.pop('net_par', None)
        selected_parasite = st.selectbox('Select a parasite to visualize the predicted PPI', parasite_list, key="net_par")

    # Set info message on initial site load
    if selected_parasite == "<select>":
        st.text('Choose 1 parasite to visualize the predicted PPI network')
    else:
        df_select = get_parasite_tissues(data_dir, selected_parasite, selected_taxids)
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
        G = generate_graph(df_select, score, web_utils.load_protein_annotations(data_dir))
            
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
        set_label_font(net)
        
        #net.show_buttons(filter_=['nodes'])
        
        
with col3:
    st.write('')


with st.container():
    if net is not None:
        st.caption('Predicted interactions between the proteins of the parasite and those of the '
                   'host, above the chosen confidence score. Each node is a protein: diamonds are '
                   'parasite proteins and circles host proteins, coloured by the organism they '
                   'belong to and drawn the larger the more central they are to the network. Each '
                   'edge is a predicted interaction, drawn the thicker the higher its confidence '
                   'score. Hover over a node for the full protein name and its identifiers.')
        html_data = ""
        # Save and read graph as HTML file (on Streamlit Sharing)
        net.save_graph(f'{path}/{selected_parasite}.html')
        with open(f'{path}/{selected_parasite}.html','r',encoding='utf-8') as HtmlFile:
            html_data = HtmlFile.read()
        # Load HTML into HTML component for display on Streamlit
        st.iframe(html_data, height=1050)
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
                            height=350
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
        enrichment = get_enrichment(df_select[df_select['weight'] >= score], data_dir)
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
                                fit_columns_on_grid_load=False,
                                enable_enterprise_modules=True,
                                height=350
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
            selected_terms = selected_rows['go_term'].values.tolist()
            enrichment_viz = enrichment_viz[enrichment_viz['go_term'].isin(selected_terms)]

        with go1:
            fig = px.scatter(enrichment_viz, x='fdr_bh', y='odds_ratio', 
                size='odds_ratio', color='go_term', height=450, 
                labels = {'fdr_bh':'FDR BH', 'odds_ratio': 'Odds ratio'})
            fig.update_traces(showlegend=False)
            st.subheader("Enriched Biological Processes -- Odds ratio vs FDR")
            st.caption('Biological processes over-represented among the host proteins of the '
                       'network. A process is the more enriched the higher up it sits (odds '
                       'ratio) and the more significant the further left (FDR).')
            st.plotly_chart(fig, width='stretch', height=400)
        with go2:
            if len(selected_terms) > 0:
                if enrichment is not None:
                    highlighted_nodes = enrichment[enrichment['go_term'].isin(selected_terms)]['nodes'].values
                    highlighted_nodes = utils.merge_list_of_lists([i.split(',') for i in highlighted_nodes])
                    highlight_color = {i: '#e7298a' for i in highlighted_nodes}
                    G = generate_graph(df_select, score, web_utils.load_protein_annotations(data_dir))
                    nx.set_node_attributes(G, "#ddd", 'color')
                    nx.set_node_attributes(G, highlight_color, 'color')
                    # Initiate PyVis network object
                    net = Network(height="450px", width="100%", bgcolor='white', font_color='#555555')
                    # Take Networkx graph and translate it to a PyVis graph format
                    net.from_nx(G)
                    G = None
                    set_label_font(net)
                    net.save_graph(f'{path}/{selected_parasite}2.html')
                    net = None
                    st.subheader("Highlighted Nodes for Selected Biological Processes")
                    st.caption('The same network, with the proteins annotated to the biological '
                               'processes selected in the table above in pink and the rest in grey.')
                    with open(f'{path}/{selected_parasite}2.html','r',encoding='utf-8') as HtmlFile:
                        html_data = HtmlFile.read()
                    st.iframe(html_data, height=500)
                    st.download_button(
                        label="Download Network as Html",
                        data=html_data,
                        file_name=f'{selected_parasite}_enrichment_network.html',
                        mime='text/html',
                    )
        
        fig = get_enrichment_summary(enrichment_table, load_ontology(data_dir))
        st.subheader("Visual Summary of Enriched Hierarchy of Biological Processes")
        st.caption('The enriched processes arranged by the Gene Ontology hierarchy, each nested '
                   'in its parent term and sized by its odds ratio, so related processes are read '
                   'together rather than as a list. Click a block to zoom in.')
        st.plotly_chart(fig, width='stretch')

st.markdown("---")
st.markdown("---")

# Footer
with st.container():
    web_utils.footer()