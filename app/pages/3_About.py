import web_utils
import streamlit as st
from css import style

style.load_css()
web_utils.show_header('About')

# Footer
with st.container():
    web_utils.footer()