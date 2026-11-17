import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import web_utils
import streamlit as st
import body_figure
import section_nav
from multi_host import (ABSENT_GROUP, NOT_EXPRESSED, NOT_TRANSFERRED, OUT_OF_REACH,
                        explain_host_specific, generate_combination_bars,
                        generate_link_matrix, get_host_coverage, get_link_combinations,
                        get_multi_host_predictions)

config = web_utils.load_config(web_utils.get_config_file())
data_dir = web_utils.get_data_dir()



# the sections on the left, with the slider that sets them under them; the title heads
# the column of figures
nav = section_nav.SectionNav('hosts_of_parasite', [
    ('Predicted interactions relative to the available host proteins', 'coverage'),
    ('Predicted interactions per infected organ', body_figure.ORGANS_ANCHOR),
    ('Shared and host-specific interactions', 'sets'),
    ('Interactions per parasite and host orthologous group', 'links')])
with nav.body:
    _, parasite = web_utils.require_selection('Hosts of a parasite', parasite=True,
                                              back=nav.back)
    st.caption('Compare the hosts the selected parasite is predicted against.')
with nav.side:
    score = web_utils.score_slider(
        help='Interactions predicted below this confidence are left out, as '
             'on the other pages. These are the smallest interactomes of the '
             'set: raising it can leave a parasite with only one host.')

with nav.body:
    df_pred = get_multi_host_predictions(data_dir, config, score)

    if parasite not in set(df_pred['taxid1_label']):
        st.info(f'{parasite} is not predicted against more than one host at this confidence.')
    else:
        edges = df_pred[df_pred['taxid1_label'] == parasite]
        all_hosts = sorted(edges['host'].unique())
        hosts_taxids = {host: tuple(sorted(set(rows['taxid2'])))
                        for host, rows in edges.groupby('host')}
        links = get_link_combinations(data_dir, config, score)[parasite]

        # the comparison first, since a host carrying more interactions may only have been
        # annotated more deeply
        nav.subheader(
            'Predicted interactions relative to the available host proteins', 'coverage',
            help='Interactions predicted in each host beside the candidate host proteins '
                 'they were drawn from: those passing the tissue and DeepLoc filters, '
                 'lying where this parasite\'s niche can reach them — the surface of the '
                 'cell for a parasite that stays outside it, the cytosol and nucleus as '
                 'well for one that gets in — and annotated to a tissue this parasite '
                 'infects. A host with more predicted interactions than another may '
                 'simply have more candidates.')
        st.caption('The interactions predicted in each host, beside how many candidate '
                   'host proteins it had.')
        st.dataframe(get_host_coverage(data_dir, config, score, parasite, hosts_taxids),
                     width='stretch', hide_index=True,
                     column_config={'candidate host proteins': st.column_config.NumberColumn(
                         help='Host proteins that pass the tissue and DeepLoc filters, '
                              'lie where this parasite can reach them, and are expressed '
                              'in a tissue it infects: the pool the predictions were '
                              'drawn from.')})

        # the figure is drawn narrower than the page, where the body outlines are read best
        nav.add('Predicted interactions per infected organ', body_figure.ORGANS_ANCHOR)
        nav.divider()
        with st.columns([2, 1])[0]:
            body_figure.show_body_figure(
                config, data_dir, edges[edges['weight'] >= score],
                tuple(taxid for host in all_hosts for taxid in hosts_taxids[host]),
                shared_color_scale=True, title_as_subheader=True)

        nav.subheader(
            f'Shared and host-specific interactions of *{parasite}*', 'sets',
            label='Shared and host-specific interactions',
            help=f'Interactions of {parasite} predicted in the set of hosts named by '
                 'the matrix below the bars. Interactions are counted as pairs of '
                 'orthologous groups rather than pairs of proteins, since the '
                 'orthologous proteins of two hosts are distinct proteins.')
        st.caption('How many interactions are found in each set of hosts.')
        st.plotly_chart(generate_combination_bars(links, all_hosts, config),
                        width='stretch')

        reasons = explain_host_specific(links, all_hosts, data_dir, config,
                                        edges['taxid1'].iloc[0], hosts_taxids)
        if reasons is not None and not reasons.empty:
            counts = reasons['Why'].value_counts()
            # `none` rather than `0`: a nought in bold reads as a value the page failed
            # to fill in
            absent = counts.get(ABSENT_GROUP, 0) or 'none'
            lead = (f'**Of the {len(reasons)} host–link combinations missing from a '
                    f'host, {absent} are missing because the host lacks the orthologous group.**')
            clauses = []
            if counts.get(NOT_EXPRESSED):
                clauses.append(f'in {counts[NOT_EXPRESSED]} the group is not annotated '
                               'to a tissue the parasite infects')
            if counts.get(OUT_OF_REACH):
                clauses.append(f'in {counts[OUT_OF_REACH]} it is expressed there, but '
                               'DeepLoc places it in compartments the parasite cannot '
                               'reach')
            if counts.get(NOT_TRANSFERRED):
                clauses.append(f'{counts[NOT_TRANSFERRED]} have a member that passes '
                               'every filter, yet no interaction is recorded — a '
                               'mismatch between the predictions and the data directory')
            if clauses:
                rest = '; '.join(clauses)
                lead += ' ' + rest[0].upper() + rest[1:] + '.'
            st.markdown(lead)
            with st.expander(f'The {len(reasons)} missing links, one row each'):
                st.dataframe(reasons, width='stretch', height=400, hide_index=True)
        elif reasons is None:
            st.caption('Run `python scripts/build_host_orthologs.py` to add the check '
                       'of whether a host missing an interaction has a protein of the '
                       'orthologous group at all.')

        nav.subheader(
            'Interactions per parasite and host orthologous group', 'links',
            help='One tile per transferred interaction: parasite orthologous groups on '
                 'the x axis, host orthologous groups on the y axis, coloured by the '
                 'host or host set that received the interaction. The groups are '
                 'labelled with the proteins they hold, with the '
                 'group ID and the full membership on hover, ordered so shared '
                 'interactions gather in the upper left. The links found in every host '
                 'are left out unless switched on below, since they are most of the '
                 'links and are counted in the bars above.')
        st.caption('Which hosts each interaction is found in, one tile per interaction.')
        # the shared links are most of the matrix, several hundred rows for Trichinella
        show_shared = st.toggle('Include links found in every host', value=False)
        shown = links if show_shared else links[links['n_hosts'] < len(all_hosts)]
        if shown.empty:
            st.text(f'Every link of {parasite} is found in every host')
        else:
            st.plotly_chart(generate_link_matrix(shown, all_hosts, config), width='stretch')

nav.finish()

st.markdown("---")


# Footer
with st.container():
    web_utils.footer()
