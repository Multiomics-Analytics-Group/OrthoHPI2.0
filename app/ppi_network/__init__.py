import os
import streamlit.components.v1 as components

# built from this file rather than the working directory, which differs between local runs
# and the Procfile
_FRONTEND_DIR = os.path.dirname(os.path.abspath(__file__))
_component = components.declare_component('ppi_network', path=_FRONTEND_DIR)


def ppi_network(nodes, edges, options, height=1000, key=None):
    '''Draws a vis.js network and reports back the interaction the user clicks on.'''
    return _component(nodes=nodes, edges=edges, options=options, height=height,
                      key=key, default=None)
