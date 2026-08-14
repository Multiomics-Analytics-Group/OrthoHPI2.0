import os
import streamlit.components.v1 as components

# The frontend is a plain index.html served by Streamlit itself, so there is no build
# step and nothing to install. The path is built from this file rather than from the
# working directory, which differs between running the app locally and the Procfile.
_FRONTEND_DIR = os.path.dirname(os.path.abspath(__file__))
_component = components.declare_component('ppi_network', path=_FRONTEND_DIR)


def ppi_network(nodes, edges, options, height=1000, key=None):
    '''
    Draws a vis.js network and reports back the interaction the user clicks on. The
    nodes, edges and options are the ones pyvis would have written into its HTML
    (Network.get_network_data()), so the network looks and behaves the same as the one
    pyvis draws -- the only thing gained is that a click reaches Python, which an HTML
    file embedded with st.iframe cannot do.

    :param list nodes: vis.js node dictionaries
    :param list edges: vis.js edge dictionaries; whatever extra keys they carry come
                       back with the click, which is how the clicked interaction is
                       identified
    :param dict options: vis.js options
    :param int height: height of the network in pixels
    :param str key: Streamlit widget key
    :return: {'nonce': ..., 'edge': {...}} of the last edge clicked, None before the
             first click. The nonce changes on every click, so clicking the same edge
             twice is told apart from a rerun that changed nothing.
    '''
    return _component(nodes=nodes, edges=edges, options=options, height=height,
                      key=key, default=None)
