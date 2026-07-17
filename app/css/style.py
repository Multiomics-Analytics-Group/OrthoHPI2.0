import os
import streamlit as st


def load_css():
    css_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "style.css")
    with open(css_path) as f:
        st.markdown('<style>{}</style>'.format(f.read()), unsafe_allow_html=True)