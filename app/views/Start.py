import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import web_utils
import organism_picker
import streamlit as st

config = web_utils.load_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()


def common(host):
    return organism_picker.split_label(host)[0].lower()


def card(slot, key, title, description, page, ready, locked_note):
    '''
    One analysis to open, as a button with what it shows under it; shaded, with what it
    still needs, until the selection allows it.
    '''
    state = 'open' if ready else 'locked'
    with slot.container(border=True, height='stretch', gap='small', key=f'card_{state}_{key}'):
        if st.button(title, key=f'open_{key}', disabled=not ready, type='primary',
                     width='stretch'):
            st.switch_page(page)
        st.caption(description if ready else locked_note)


web_utils.show_header()

host, parasite = web_utils.get_selection()
parasite_hosts = organism_picker.hosts_of(config, data_dir, parasite) if parasite else []

general, of_host, of_parasite, network = st.columns(4)
card(general, 'general', 'General analysis',
     'Every prediction, across all hosts and parasites.',
     web_utils.GENERAL, True, '')
card(of_host, 'of_host', 'Compare parasites of this host',
     f'Compare the parasites of {common(host)}.' if host else '',
     web_utils.PARASITES_OF_HOST, bool(host), 'Pick a host below.')
card(of_parasite, 'of_parasite', 'Compare hosts of this parasite',
     f'Compare *{parasite}* across {", ".join(common(h) for h in parasite_hosts)}.'
     if parasite else '',
     web_utils.HOSTS_OF_PARASITE, len(parasite_hosts) > 1,
     f'*{parasite}* has predictions in one host only.' if parasite
     else 'Pick a parasite found in several hosts below.')
card(network, 'network', 'View network',
     f'The predicted network of *{parasite}* in {common(host)}.' if host and parasite else '',
     web_utils.NETWORK, bool(host and parasite), 'Pick a host and a parasite below.')

if host or parasite:
    with st.container(horizontal=True, vertical_alignment='center', key='selection_line'):
        parts = [f'Host: **{host}**' if host else 'No host',
                 f'Parasite: ***{parasite}***' if parasite else 'no parasite']
        st.markdown(' · '.join(parts))
        if st.button('Clear', icon=':material/close:', type='tertiary'):
            web_utils.set_selection(None, None)
            st.rerun()

organism_picker.show(config, data_dir)

st.page_link(web_utils.ABOUT, label='About OrthoHPI', icon=':material/info:')
st.markdown('---')
with st.container():
    web_utils.footer()
