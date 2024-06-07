import utils
import web_utils
from css import style
import pandas as pd
import streamlit as st
from streamlit_extras.switch_page_button import switch_page
from stmol import showmol
import structure_visualizer as strv
from st_aggrid import GridOptionsBuilder, AgGrid

style.load_css()
page = web_utils.show_pages_menu(index=2)
if page == "Home":
    switch_page("orthohpi home")
elif page == "Predicted Host-parasite PPIs":
    switch_page("predicted host-parasite ppis")
elif page == "About":
    switch_page('about')


def get_structures(query_proteins):
    structures = strv.get_alphafold_structure(query_proteins=query_proteins)
    
    return structures

def show_structure(pdb_file):
    xyzview = strv.generate_mol_structure(pdb_file=pdb_file)
    showmol(xyzview, height = 500,width=700)


config = utils.read_config('config.yml')
predictions = utils.read_parquet_file(input_file='data/annotated_predictions.parquet')
predictions['weight'] = predictions['weight'].astype(float)
tissues = utils.read_parquet_file(input_file='data/tissues_cell_types.parquet')
pred_tissues = pd.merge(predictions, tissues.rename({'Gene': 'target'}, axis=1), on='target', how='left')
tissues = None
predictions = None
parasite_list = ['<select>'] + pred_tissues['taxid1_label'].sort_values().unique().tolist()


st.markdown("<h1 style='text-align: center; color: #023858;'>OrthoHPI 2.0</h1>", unsafe_allow_html=True)
st.markdown("<h3 style='text-align: center; color: #2b8cbe;'>AlphaFold Predicted Structure Interactors</h3>", unsafe_allow_html=True)


col1, col2, col3 = st.columns(3)

with col1:
    st.write('')


with col2:
    # Implement multiselect dropdown menu for option selection
    selected_parasite = st.selectbox('Select a parasite to visualize the predicted PPI', parasite_list, key="struct_par")
    if selected_parasite == "<select>":
        st.text('Choose 1 parasite to visualize the predicted PPI network')
        selected_cols = None
    else:
        selected_cols = ['taxid1_label', 'source_name', 'source', 'taxid2_label',
        'target_name', 'target', 'experimental_evidence_score', 'databases_evidence_score',
        'weight', 'group1', 'group2', 'source_uniprot', 'target_uniprot']
        score = st.slider('Confidence score', 0.4, 0.9, 0.7)
  
with col3:
    st.write('')
        
        
if selected_cols is not None:    
    with st.container():
        df_select = pred_tissues.loc[pred_tissues['taxid1_label'] == selected_parasite]
        df_select = web_utils.filter_tissues(config, df_select)
        df_select = df_select[df_select['weight'] >= score]
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
                            update_mode='MODEL_CHANGED',
                            fit_columns_on_grid_load=False,
                            enable_enterprise_modules=True,
                            height=350,
                            reload_data=False
                        )
        selected_rows = grid_response['selected_rows']
        if selected_rows is not None and len(selected_rows) > 0:
            query_proteins = dict(selected_rows[['source_name', 'source_uniprot']].values)
            query_proteins.update(dict(selected_rows[['target_name', 'target_uniprot']].values))
            structures = get_structures(query_proteins)
            cols = st.columns(2)
            i = 0
            for protein in structures:
                pdb_file, url, website = structures[protein]
                with cols[i]:
                    if pdb_file is not None:
                        st.markdown(f'''<h4>AlphaFold structure {protein}</h4>''',
                        unsafe_allow_html=True)
                        show_structure(pdb_file=pdb_file)
                        st.markdown(f'''
                                <a href={url}><button>PDB file</button></a>
                                <a href={website}><button>AlphaFold EBI</button></a>''',
                                unsafe_allow_html=True)
                    else:
                        st.markdown(f'''<h4>AlphaFold structure {protein}</h4>''',
                        unsafe_allow_html=True)
                        st.markdown(f'''<h5>AlphaFold prediction Not Available</h5>''',
                        unsafe_allow_html=True)
                    i += 1


st.markdown("---")
st.markdown("---")

# Footer
with st.container():
    web_utils.footer()