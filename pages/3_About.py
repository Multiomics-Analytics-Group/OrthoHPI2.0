import web_utils
import streamlit as st

st.set_page_config(layout="wide", page_title="About", menu_items={})

# Footer
with st.container():
    web_utils.footer()