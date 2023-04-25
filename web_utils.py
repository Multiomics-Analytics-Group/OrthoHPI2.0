import streamlit as st
from streamlit_option_menu import option_menu

def show_pages_menu(index=0):
    selected = option_menu(
    menu_title=None,  # required
    options=["Home", "Predicted Host-parasite PPIs", "Predicted PPI structures", "About"],  # required
    icons=["house", "diagram-3", "cast", "chat-text"],  # optional
    menu_icon="cast",  # optional
    default_index=index,  # optional
    orientation="horizontal",
    )
    return selected


def filter_tissues(config, df):
    source = df['taxid1'].unique()[0]
    mapped_tissues = config['tissues']
    tissues = [mapped_tissues[t].lower() for t in config['parasites'][int(source)]['tissues']]
    df = df[df['Tissue'].isin(tissues)]
    
    return df

def footer():
    st.write("Developed with data from:")

    cols = st.columns(6)
    with cols[0]:
        st.image('images/eggnog.png', width=200)
    with cols[1]:
        st.image('images/string.png', width=200)
    with cols[2]:
        st.image('images/hpa.png', width=200)
    with cols[3]:
        st.image('images/tissues.png', width=200)
    with cols[4]:
        st.image('images/compartments.png', width=200)
    with cols[5]:
        st.image('images/ebi.png', width=200)

    st.write("Code available at: https://github.com/Multiomics-Analytics-Group/OrthoHPI2.0")