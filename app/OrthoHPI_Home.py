import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import web_utils
import streamlit as st
from css import style

st.set_page_config(layout="wide", page_title="OrthoHPI 2.0", menu_items={})
st.session_state.data_dir = 'data'
st.session_state.config_file = 'config.yml'
style.load_css()
st.logo(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logo.svg'), size='large')

page = web_utils.navigation()
page.run()
