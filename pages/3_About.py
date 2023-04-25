import web_utils
import streamlit as st
from streamlit_extras.switch_page_button import switch_page
from css import style

style.load_css()
page = web_utils.show_pages_menu(index=3)
if page == "Home":
    switch_page("orthohpi home")
elif page == "Predicted Host-parasite PPIs":
    switch_page("predicted host-parasite ppis")
elif page == "Predicted PPI structures":
    switch_page('interaction structures')

# Footer
with st.container():
    web_utils.footer()