import streamlit as st

def footer():
    st.write("Developed with data from:")

    cols = st.columns(5)
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

    st.write("Code available at: https://github.com/Multiomics-Analytics-Group/OrthoHPI2.0")