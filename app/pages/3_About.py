import web_utils
import streamlit as st
from css import style

st.set_page_config(layout="wide", page_title="OrthoHPI 2.0", menu_items={})
style.load_css()
web_utils.show_header('About')

# Footer
with st.container():
    web_utils.footer()