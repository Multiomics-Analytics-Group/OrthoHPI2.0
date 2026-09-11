"""
Draws the hosts and parasites of the study on the home page: a card per host with the
number of parasites it has predictions for, and a panel per parasite group listing its
species, each marked with a dot per host it infects.

Everything is read from the configuration, so the overview follows the study whenever a
host or a parasite is added. It is the same picture as docs/hosts_parasites.svg, drawn in
HTML rather than as a fixed 1600 x 900 frame so that it wraps to the page and sits among
the other figures; the group order, the plain-language names of the groups and their icons
are imported from the script that writes the slide, so the two cannot drift apart.
"""
import html
import re

import streamlit as st

from scripts.build_hosts_parasites_figure import GROUP_ICON, GROUP_ORDER, GROUP_SUB, ICONS

# `Homo sapiens (human)`: the latin name and, in parentheses, the one the cards are headed with
LABEL = re.compile(r'^(?P<latin>[^(]+?)\s*\((?P<common>[^)]+)\)\s*$')

# the panel Nematoda takes up: twice the width and height of the others, with its names in
# two columns, since it holds more species than the next three groups together
WIDE_GROUP = 'Nematoda'

STYLE = '''
<style>
.hp { font-size: 0.95rem; }
.hp-label { font-size: 0.75rem; font-weight: 600; letter-spacing: 0.1em; opacity: 0.6;
            margin: 1.2rem 0 0.5rem; }
.hp-hosts { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 1rem; }
.hp-host { position: relative; border-radius: 14px; padding: 0.9rem 1rem 0.9rem 1.5rem;
           overflow: hidden; }
.hp-host::before { content: ""; position: absolute; inset: 0; background: var(--c);
                   opacity: 0.08; }
.hp-host::after { content: ""; position: absolute; left: 0; top: 14px; bottom: 14px; width: 5px;
                  border-radius: 2.5px; background: var(--c); }
.hp-host b { font-size: 1.3rem; font-weight: 600; display: block; }
.hp-host i { display: block; opacity: 0.65; font-size: 0.9rem; }
.hp-host span { display: block; color: var(--c); font-weight: 600; margin-top: 0.5rem; }
.hp-groups { display: grid; grid-template-columns: repeat(auto-fit, minmax(230px, 1fr));
             gap: 1rem; align-items: start; }
.hp-group { position: relative; border-radius: 14px; padding: 1rem 0.75rem 0.7rem;
            background: color-mix(in srgb, currentColor 4%, transparent); overflow: hidden; }
.hp-group::before { content: ""; position: absolute; left: 0; right: 0; top: 0; height: 5px;
                    background: var(--c); }
.hp-group.wide { grid-column: span 2; }
.hp-group.wide ul { columns: 2; column-gap: 1rem; }
/* on a page five panels wide the seven fill two rows exactly if Nematoda takes both;
   with fewer across, the span leaves a hole under it instead */
@media (min-width: 1400px) { .hp-group.wide { grid-row: span 2; } }
@media (max-width: 560px) { .hp-group.wide { grid-column: auto; }
                            .hp-group.wide ul { columns: 1; } }
.hp-head { display: flex; align-items: center; gap: 0.7rem; margin-bottom: 0.6rem; }
.hp-head svg { flex: none; width: 40px; height: 40px; }
.hp-head b { display: block; font-size: 1.1rem; font-weight: 600; }
.hp-head small { display: block; opacity: 0.65; font-size: 0.8rem; }
.hp-head em { margin-left: auto; font-style: normal; font-size: 1.35rem; font-weight: 700;
              color: var(--c); }
.hp ul { list-style: none; margin: 0; padding: 0; }
.hp li { display: flex; align-items: flex-start; gap: 0.1rem; padding: 0.3rem 0;
         break-inside: avoid; }
.hp li i { margin-left: 0.35rem; font-size: 0.85rem; line-height: 1.3; min-width: 0; }
.hp-dot { flex: none; width: 9px; height: 9px; border-radius: 50%; background: var(--c); }
/* a name that wraps keeps its dots on the first line, centred on it */
.hp li .hp-dot { margin-top: calc(0.85rem * 0.65 - 4.5px); }
.hp-legend { display: flex; flex-wrap: wrap; align-items: center; gap: 0.3rem 0.9rem;
             margin-top: 0.9rem; font-size: 0.8rem; opacity: 0.75; }
.hp-legend span { display: inline-flex; align-items: center; gap: 0.35rem; }
</style>
'''


def split_label(label):
    '''
    The common and the latin name of a host, from the label the configuration gives it.

    :param str label: e.g. `Homo sapiens (human)`
    :return: (`Human`, `Homo sapiens`); a label without parentheses is both
    '''
    found = LABEL.match(label)
    if not found:
        return label, label

    return found['common'].capitalize(), found['latin']


def dots(host_ids, hosts, infected):
    '''
    One slot per host, in the order of the configuration, so that the dots of every
    species line up and the position of a dot says which host it is before its colour does.

    :param list host_ids: the hosts of the study
    :param dict hosts: their configuration
    :param set infected: the hosts of the species
    :return: the html of the slots
    '''
    return ''.join(
        f'<span class="hp-dot" style="--c: {hosts[h]["color"]}" title="{split_label(hosts[h]["label"])[0]}"></span>'
        if h in infected else '<span class="hp-dot" style="visibility: hidden"></span>'
        for h in host_ids)


def host_card(taxid, host, count):
    common, latin = split_label(host['label'])
    return (f'<div class="hp-host" style="--c: {host["color"]}"><b>{html.escape(common)}</b>'
            f'<i>{html.escape(latin)}</i><span>{count} parasite{"s" * (count != 1)}</span></div>')


def group_panel(group, color, members, host_ids, hosts):
    icon = f'<svg viewBox="0 0 100 100">{ICONS[GROUP_ICON[group]](color)}</svg>'
    items = ''.join(
        f'<li>{dots(host_ids, hosts, infected)}<i>{html.escape(name)}</i></li>'
        for name, infected in members)
    wide = ' wide' if group == WIDE_GROUP else ''
    return (f'<div class="hp-group{wide}" style="--c: {color}"><div class="hp-head">{icon}'
            f'<div><b>{html.escape(group)}</b><small>{html.escape(GROUP_SUB.get(group, ""))}</small></div>'
            f'<em>{len(members)}</em></div><ul>{items}</ul></div>')


def show(config, interactions=None):
    '''
    Draws the overview: the hosts across the top, the parasites grouped below them, and a
    legend for the dots.

    :param dict config: the configuration of the study
    :param int interactions: number of predicted interactions, counted once per host a
                             parasite infects, to state in the caption with the species
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
    panels = ''.join(
        group_panel(group, palette.get(group, '#999999'), sorted(by_group[group]), host_ids, hosts)
        for group in groups)
    legend = ''.join(
        f'<span><span class="hp-dot" style="--c: {host["color"]}"></span>'
        f'{html.escape(split_label(host["label"])[0])}</span>'
        for host in hosts.values())

    caption = f'{len(parasites)} parasite species in {len(hosts)} host species'
    if interactions is not None:
        caption += f', {interactions:,} predicted host-parasite interactions'
    st.caption(caption + '.')
    st.markdown(
        STYLE + f'<div class="hp"><div class="hp-label">HOSTS</div><div class="hp-hosts">{cards}</div>'
        f'<div class="hp-label">PARASITES</div><div class="hp-groups">{panels}</div>'
        f'<div class="hp-legend"><span>Host:</span>{legend}</div></div>',
        unsafe_allow_html=True)
