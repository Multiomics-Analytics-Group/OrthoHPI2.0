import sys, os, base64
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import web_utils
import streamlit as st
from css import style

APP_DIR = os.path.dirname(os.path.abspath(__file__))
FAVICON = os.path.join(APP_DIR, "favicon.png")
st.set_page_config(layout="wide", page_title="OrthoHPI 2.0", menu_items={}, page_icon=FAVICON)
st.session_state.data_dir = 'data'
st.session_state.config_file = 'config.yml'
style.load_css()
st.logo(os.path.join(APP_DIR, "logo.svg"), size="large")

# Streamlit swaps the address of the page's icon link, which Firefox ignores; a new link it
# does pick up
with open(FAVICON, 'rb') as f:
    favicon_uri = 'data:image/png;base64,' + base64.b64encode(f.read()).decode()
st.html(f"""<script>
if (!document.getElementById('orthohpi-icon')) {{
    document.querySelectorAll("link[rel~='icon']").forEach(link => link.remove());
    const link = document.createElement('link');
    link.id = 'orthohpi-icon';
    link.rel = 'icon';
    link.type = 'image/png';
    link.href = '{favicon_uri}';
    document.head.appendChild(link);
}}
</script>""", unsafe_allow_javascript=True)

page = web_utils.navigation()
page.run()
