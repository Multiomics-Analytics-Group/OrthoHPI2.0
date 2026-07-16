import web_utils
import streamlit as st
from css import style

style.load_css()
page = web_utils.show_pages_menu(index=3)
if page == "Home":
    st.switch_page("OrthoHPI_Home.py")
elif page == "Predicted Host-parasite PPIs":
    st.switch_page("pages/1_Predicted_Host-Parasite_PPIs.py")
elif page == "Predicted PPI structures":
    st.switch_page('pages/2_Interaction_structures.py')

# Footer
with st.container():
    web_utils.footer()