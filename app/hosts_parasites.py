'''
Draws the hosts and parasites of the study on the home page: a chip per host with the number
of parasites it has predictions for, and a box per parasite group listing its species,
each marked with a dot per host it infects.
'''
import html
import re

import streamlit as st

from scripts.build_hosts_parasites_figure import GROUP_ICON, GROUP_ORDER, GROUP_SUB, ICONS

# `Homo sapiens (human)`: the latin name and, in parentheses, the one the cards are headed
# with
LABEL = re.compile(r'^(?P<latin>[^(]+?)\s*\((?P<common>[^)]+)\)\s*$')

# Nematoda takes the wide box, with its names in two columns
WIDE_GROUP = 'Nematoda'

STYLE = '''
<style>
.hp { font-size: 0.95rem; }
.hp-label { font-size: 0.75rem; font-weight: 600; letter-spacing: 0.1em; opacity: 0.6;
            margin: 1.2rem 0 0.5rem; }
.hp-hosts { display: flex; flex-wrap: wrap; gap: 0.5rem; }
.hp-host { display: inline-flex; align-items: baseline; gap: 0.45rem; border-radius: 999px;
           padding: 0.3rem 0.9rem; background: color-mix(in srgb, var(--c) 10%, transparent); }
.hp-host .hp-dot { align-self: center; }
.hp-host b { font-weight: 600; }
.hp-host i { opacity: 0.65; font-size: 0.85rem; }
.hp-host span { color: var(--c); font-weight: 700; }
/* Nematoda, the largest group, takes a box with its names in two columns; the other groups
   pack into columns beside it, or below it on a narrow page */
.hp-groups { display: flex; flex-wrap: wrap; gap: 0.7rem; align-items: flex-start; }
.hp-rest { flex: 3 1 460px; columns: 210px; column-gap: 0.7rem; }
.hp-group { position: relative; border-radius: 10px; padding: 0.65rem 0.7rem 0.45rem;
            background: color-mix(in srgb, currentColor 4%, transparent); overflow: hidden;
            break-inside: avoid; margin-bottom: 0.7rem; }
.hp-group::before { content: ""; position: absolute; left: 0; right: 0; top: 0; height: 4px;
                    background: var(--c); }
.hp-group.wide { flex: 2 1 420px; margin-bottom: 0; }
.hp-group.wide ul { columns: 2; column-gap: 1rem; }
@media (max-width: 560px) { .hp-group.wide ul { columns: 1; } }
.hp-head { display: flex; align-items: center; gap: 0.1rem 0.45rem; margin-bottom: 0.3rem; }
/* the subtitle goes under the name when the box is too narrow for both */
.hp-head div { display: flex; flex-wrap: wrap; align-items: baseline; column-gap: 0.45rem; }
.hp-head svg { flex: none; width: 22px; height: 22px; }
.hp-head b { font-weight: 600; }
.hp-head small { opacity: 0.65; font-size: 0.75rem; white-space: nowrap; }
.hp-head em { margin-left: auto; font-style: normal; font-weight: 700; color: var(--c); }
.hp ul { list-style: none; margin: 0; padding: 0; }
.hp li { display: flex; align-items: center; gap: 0.15rem; padding: 0.1rem 0;
         break-inside: avoid; }
.hp li i { margin-left: 0.3rem; font-size: 0.84rem; }
.hp-dot { flex: none; width: 8px; height: 8px; border-radius: 50%; background: var(--c); }
</style>
'''


def split_label(label):
    '''
    The common and the latin name of a host, from the label the configuration gives it.
    '''
    found = LABEL.match(label)
    if not found:
        return label, label

    return found['common'].capitalize(), found['latin']


def dots(host_ids, hosts, infected):
    '''
    A dot per host the species infects, in the order of the configuration.
    '''
    return ''.join(
        f'<span class="hp-dot" style="--c: {hosts[h]["color"]}" title="{split_label(hosts[h]["label"])[0]}"></span>'
        for h in host_ids if h in infected)


def host_card(taxid, host, count):
    common, latin = split_label(host['label'])
    return (f'<div class="hp-host" style="--c: {host["color"]}"><span class="hp-dot"></span>'
            f'<b>{html.escape(common)}</b><i>{html.escape(latin)}</i>'
            f'<span title="parasites">{count}</span></div>')


def group_panel(group, color, members, host_ids, hosts):
    icon = f'<svg viewBox="0 0 100 100">{ICONS[GROUP_ICON[group]](color)}</svg>'
    items = ''.join(
        f'<li>{dots(host_ids, hosts, infected)}<i>{html.escape(name)}</i></li>'
        for name, infected in members)
    wide = ' wide' if group == WIDE_GROUP else ''
    return (f'<div class="hp-group{wide}" style="--c: {color}"><div class="hp-head">{icon}'
            f'<div><b>{html.escape(group)}</b><small>{html.escape(GROUP_SUB.get(group, ""))}</small></div>'
            f'<em>{len(members)}</em></div><ul>{items}</ul></div>')


def show(config):
    '''
    Draws the overview: the hosts as a row of chips, which also key the dots, and a box
    per parasite group listing its species.
    '''
    hosts = config['hosts']
    parasites = config['parasites']
    palette = config.get('parasite_groups', {})
    host_ids = list(hosts)

    by_group = {}
    for parasite in parasites.values():
        by_group.setdefault(parasite['group'], []).append(
            (parasite['label'], set(parasite.get('hosts', []))))
    # the groups the slide orders, then any the configuration has since gained
    groups = [g for g in GROUP_ORDER if g in by_group] + sorted(set(by_group) - set(GROUP_ORDER))

    cards = ''.join(
        host_card(taxid, host, sum(1 for p in parasites.values() if taxid in p.get('hosts', [])))
        for taxid, host in hosts.items())
    panels = {
        group: group_panel(group, palette.get(group, '#999999'), sorted(by_group[group]),
                           host_ids, hosts)
        for group in groups}
    rest = ''.join(panel for group, panel in panels.items() if group != WIDE_GROUP)
    st.markdown(
        STYLE + f'<div class="hp"><div class="hp-label">HOSTS</div><div class="hp-hosts">{cards}</div>'
        f'<div class="hp-label">PARASITES</div><div class="hp-groups">'
        f'{panels.get(WIDE_GROUP, "")}<div class="hp-rest">{rest}</div></div></div>',
        unsafe_allow_html=True)
