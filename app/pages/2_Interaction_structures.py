import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import utils
import web_utils
from css import style
import pandas as pd
import streamlit as st
import structure_visualizer as strv
from st_aggrid import GridOptionsBuilder, AgGrid

style.load_css()
page = web_utils.show_pages_menu(index=2)
if page == "Home":
    st.switch_page("OrthoHPI_Home.py")
elif page == "Predicted Host-parasite PPIs":
    st.switch_page("pages/1_Predicted_Host-Parasite_PPIs.py")
elif page == "About":
    st.switch_page('pages/3_About.py')


def get_structures(query_proteins):
    structures = strv.get_alphafold_structure(query_proteins=query_proteins)
    
    return structures

def show_structure(pdb_file):
    xyzview = strv.generate_mol_structure(pdb_file=pdb_file)
    # same as stmol.showmol, which still embeds through the deprecated st.components.v1.html
    st.iframe(xyzview._make_html(), height=500, width=700)


config = utils.read_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()
predictions = utils.read_parquet_file(input_file=f'{data_dir}/predictions.parquet')
predictions['weight'] = predictions['weight'].astype(float)
tissues = utils.read_parquet_file(input_file=f'{data_dir}/tissues_cell_types.parquet')
pred_tissues = pd.merge(predictions, tissues.rename({'Gene': 'target'}, axis=1), on='target', how='left')
tissues = None
predictions = None


st.markdown("<h1 style='text-align: center; color: #023858;'>OrthoHPI 2.0</h1>", unsafe_allow_html=True)
st.markdown("<h3 style='text-align: center; color: #2b8cbe;'>AlphaFold Predicted Structure Interactors</h3>", unsafe_allow_html=True)


col1, col2, col3 = st.columns(3)

with col1:
    st.write('')


with col2:
    # the host carries over from whichever page it was last chosen on
    selected_host, selected_taxids = web_utils.host_selector(config, pred_tissues, 'Select a host')

    if selected_host == web_utils.NO_HOST:
        st.text('Choose 1 host to explore the predicted host-parasite interactions')
        selected_parasite = "<select>"
    else:
        # only parasites that infect the selected host
        host_pred = pred_tissues[pred_tissues['taxid2'].isin(selected_taxids)]
        parasite_list = ['<select>'] + host_pred['taxid1_label'].sort_values().unique().tolist()
        # switching host can leave a parasite selected that the new host does not have
        if st.session_state.get('struct_par') not in parasite_list:
            st.session_state.pop('struct_par', None)
        selected_parasite = st.selectbox('Select a parasite to visualize the predicted PPI', parasite_list, key="struct_par")

    if selected_parasite == "<select>":
        st.text('Choose 1 parasite to visualize the predicted PPI network')
        selected_cols = None
    else:
        selected_cols = ['taxid1_label', 'source_name', 'source_full_name', 'source',
        'taxid2_label', 'target_name', 'target_full_name', 'target',
        'experimental_evidence_score', 'databases_evidence_score',
        'weight', 'group1', 'group2', 'source_uniprot', 'target_uniprot']
        score = st.slider('Confidence score', 0.4, 0.9, 0.7)
  
with col3:
    st.write('')
        
        
if selected_cols is not None:    
    with st.container():
        df_select = host_pred.loc[host_pred['taxid1_label'] == selected_parasite]
        df_select = web_utils.filter_tissues(config, df_select)
        df_select = df_select[df_select['weight'] >= score]
        # the descriptive protein name STRING carries, which the predictions only keep
        # the short version of; empty for the proteins the annotations do not cover
        annotations = web_utils.load_protein_annotations(data_dir)
        df_select = df_select.assign(
            source_full_name=df_select['source'].map(annotations).fillna(''),
            target_full_name=df_select['target'].map(annotations).fillna(''))
        df_select = df_select[selected_cols].drop_duplicates(['source_name', 'target_name'])


        gb = GridOptionsBuilder.from_dataframe(df_select)
        gb.configure_pagination(paginationAutoPageSize=True) #Add pagination
        gb.configure_side_bar() #Add a sidebar
        gb.configure_selection('single', use_checkbox=True, 
                groupSelectsChildren="Group checkbox select children") #Enable multi-row selection
        gridOptions = gb.build()
        grid_response = AgGrid(
                            df_select,
                            gridOptions=gridOptions,
                            data_return_mode='AS_INPUT',
                            fit_columns_on_grid_load=False,
                            enable_enterprise_modules=True,
                            height=350
                        )
        selected_rows = grid_response['selected_rows']
        if selected_rows is not None and len(selected_rows) > 0:
            query_proteins = dict(selected_rows[['source_name', 'source_uniprot']].values)
            query_proteins.update(dict(selected_rows[['target_name', 'target_uniprot']].values))
            structures = get_structures(query_proteins)
            st.caption('AlphaFold model of each of the two proteins of the interaction selected '
                       'above, shown as a cartoon of the backbone and coloured from the N to the '
                       'C terminus. The two are predicted separately, so the models show the '
                       'shape of each partner, not how they dock onto each other.')
            cols = st.columns(2)
            for i, protein in enumerate(structures):
                pdb_file, url, website = structures[protein]
                with cols[i % len(cols)]:
                    st.markdown(f'''<h4>AlphaFold structure {protein}</h4>''',
                    unsafe_allow_html=True)
                    if pdb_file is not None:
                        show_structure(pdb_file=pdb_file)
                        st.markdown(f'''
                                <a href={url}><button>PDB file</button></a>
                                <a href={website}><button>AlphaFold EBI</button></a>''',
                                unsafe_allow_html=True)
                    else:
                        st.markdown(f'''<h5>AlphaFold prediction Not Available</h5>''',
                        unsafe_allow_html=True)


st.markdown("---")
st.markdown("---")

# Footer
with st.container():
    web_utils.footer()