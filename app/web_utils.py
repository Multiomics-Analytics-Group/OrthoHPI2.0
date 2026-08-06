import streamlit as st
from streamlit_option_menu import option_menu

def show_pages_menu(index=0):
    selected = option_menu(
    menu_title=None,  # required
    options=["Home", "Predicted Host-parasite PPIs", "Predicted PPI structures", "About"],  # required
    icons=["house", "diagram-3", "cast", "chat-text"],  # optional
    menu_icon="cast",  # optional
    default_index=index,  # optional
    orientation="horizontal",
    )
    return selected

def get_data_dir():
    return st.session_state.get('data_dir', 'data')

def get_config_file():
    return st.session_state.get('config_file', 'config.yml')


def get_host_groups(config, predictions):
    '''
    Maps each host group offered in the app to the taxids it covers. Hosts sharing a
    `group` in the config (rat and mouse -> Rodent) are presented as one option, hosts
    without one stand alone under their label. Groups with no predicted interaction are
    dropped so the selector never offers an empty host.

    :param dict config: parsed configuration
    :param predictions: predictions dataframe, used to drop groups without predictions
    :return: {group label: [taxid as str, ...]} sorted by group label
    '''
    predicted = set(predictions['taxid2'])
    groups = {}
    for taxid, host in config['hosts'].items():
        taxid = str(taxid)
        if taxid in predicted:
            groups.setdefault(host.get('group', host['label']), []).append(taxid)

    return {group: groups[group] for group in sorted(groups)}


# The host chosen on any page applies to all of them. Streamlit keeps session_state
# across pages, but widget state is tied to the widget, so the choice is kept in a
# key of our own and each page's selectbox is seeded from it and writes back. That
# also means the selectboxes need no key and never collide.
HOST_STATE_KEY = 'selected_host'
NO_HOST = '<select>'


def host_selector(config, predictions, label='Select a host'):
    '''
    Draws the host selectbox, shared across pages through HOST_STATE_KEY.

    :param dict config: parsed configuration
    :param predictions: predictions dataframe, used to build the options
    :param str label: label shown above the selectbox
    :return: (group label, tuple of taxids as str); NO_HOST and () when nothing is chosen
    '''
    groups = get_host_groups(config, predictions)
    options = [NO_HOST] + list(groups)

    current = st.session_state.get(HOST_STATE_KEY, NO_HOST)
    if current not in options:
        # the stored host has no predictions in this data_dir (the snapshot entrypoints
        # each point at their own), so fall back rather than raise on index()
        current = NO_HOST

    selected = st.selectbox(label, options, index=options.index(current))
    st.session_state[HOST_STATE_KEY] = selected

    return selected, tuple(groups.get(selected, ()))


def filter_tissues(config, df):
    source = df['taxid1'].unique()[0]
    mapped_tissues = config['tissues']
    tissues = [mapped_tissues[t].lower() for t in config['parasites'][int(source)]['tissues']]
    df = df[df['Tissue'].isin(tissues)]
    
    return df

def footer():
    st.write("Developed with data from:")

    cols = st.columns(6)
    with cols[0]:
        st.image('images/eggnog.png', width=200)
    with cols[1]:
        st.image('images/string.png', width=200)
    with cols[2]:
        st.image('images/hpa.png', width=200)
    with cols[3]:
        st.image('images/tissues.png', width=200)
    with cols[4]:
        st.markdown("[DeepLoc 2.0](https://services.healthtech.dtu.dk/services/DeepLoc-2.0/)")
    with cols[5]:
        st.image('images/ebi.png', width=200)

    st.write("Code available at: https://github.com/Multiomics-Analytics-Group/OrthoHPI2.0")