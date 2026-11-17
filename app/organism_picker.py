'''
The hosts and parasites of the study, drawn on the start page as the way to pick them: a
chip per host with the number of parasites it has predictions for, and a box per parasite
group listing its species, each marked with a dot per host it infects. Hovering a host
lights up its parasites, hovering a parasite its hosts, and a click selects either.
'''
import html
import re

import streamlit as st

import web_utils
from scripts.build_hosts_parasites_figure import GROUP_ICON, GROUP_ORDER, GROUP_SUB, ICONS

# `Homo sapiens (human)`: the latin name and, in parentheses, the one the cards are headed
# with
LABEL = re.compile(r'^(?P<latin>[^(]+?)\s*\((?P<common>[^)]+)\)\s*$')

# Nematoda takes the wide box, a row of its own with its names in columns
WIDE_GROUP = 'Nematoda'
# separates the hosts of a parasite in its data-hosts attribute
SEP = '|'
# session key of the component, read by its click callback
PICKER_KEY = 'organism_picker'

CSS = '''
.hp { font-size: 0.95rem; }
.hp-label { font-size: 0.75rem; font-weight: 600; letter-spacing: 0.1em; opacity: 0.6;
            margin: 1.2rem 0 0.5rem; }
.hp-label:first-child { margin-top: 0; }
.hp-hosts { display: flex; flex-wrap: wrap; gap: 0.5rem; }
.hp-host { display: inline-flex; align-items: baseline; gap: 0.45rem; border-radius: 999px;
           padding: 0.3rem 0.9rem; background: color-mix(in srgb, var(--c) 10%, transparent);
           border: 2px solid transparent; cursor: pointer; user-select: none;
           transition: opacity 0.15s, background 0.15s, border-color 0.15s; }
.hp-host .hp-dot { align-self: center; }
.hp-host b { font-weight: 600; }
.hp-host i { opacity: 0.65; font-size: 0.85rem; }
.hp-host span { color: var(--c); font-weight: 700; }
.hp-host:hover { background: color-mix(in srgb, var(--c) 20%, transparent); }
.hp-host.sel { border-color: var(--c); background: color-mix(in srgb, var(--c) 22%, transparent); }
/* Nematoda, the largest group, takes a row of its own with its names in as many columns as
   fit; the other groups pack into columns below it */
.hp-groups { display: flex; flex-direction: column; gap: 0.7rem; }
.hp-rest { columns: 210px; column-gap: 0.7rem; }
.hp-group { position: relative; border-radius: 10px; padding: 0.65rem 0.7rem 0.45rem;
            background: color-mix(in srgb, currentColor 4%, transparent); overflow: hidden;
            break-inside: avoid; margin-bottom: 0.7rem; }
.hp-group::before { content: ""; position: absolute; left: 0; right: 0; top: 0; height: 4px;
                    background: var(--c); }
.hp-group.wide { margin-bottom: 0; }
.hp-group.wide ul { columns: 210px; column-gap: 1rem; }
.hp-head { display: flex; align-items: center; gap: 0.1rem 0.45rem; margin-bottom: 0.3rem; }
/* the subtitle goes under the name when the box is too narrow for both */
.hp-head div { display: flex; flex-wrap: wrap; align-items: baseline; column-gap: 0.45rem; }
.hp-head svg { flex: none; width: 22px; height: 22px; }
.hp-head b { font-weight: 600; }
.hp-head small { opacity: 0.65; font-size: 0.75rem; white-space: nowrap; }
.hp-head em { margin-left: auto; font-style: normal; font-weight: 700; color: var(--c); }
.hp ul { list-style: none; margin: 0; padding: 0; }
.hp li { display: flex; align-items: center; gap: 0.15rem; padding: 0.1rem 0.3rem;
         margin: 0 -0.3rem; border-radius: 6px; break-inside: avoid; cursor: pointer;
         transition: opacity 0.15s, background 0.15s; }
.hp li i { margin-left: 0.3rem; font-size: 0.84rem; }
.hp li:hover { background: color-mix(in srgb, currentColor 8%, transparent); }
.hp li.sel { background: color-mix(in srgb, currentColor 14%, transparent); font-weight: 600; }
.hp-dot { flex: none; width: 8px; height: 8px; border-radius: 50%; background: var(--c); }
/* a parasite without predictions, or outside the selected host, cannot be picked */
.hp li.na, .hp li.off { opacity: 0.3; cursor: default; }
.hp li.na:hover, .hp li.off:hover { background: none; }
/* while a host is hovered, the parasites it does not have fade; while a parasite is, the
   hosts it does not infect */
.hp.hovering-host .hp-host:not(.lit), .hp.hovering-host li:not(.lit),
.hp.hovering-parasite .hp-host:not(.lit) { opacity: 0.25; }
.hp.hovering-host li.lit, .hp.hovering-host .hp-host.lit,
.hp.hovering-parasite .hp-host.lit { opacity: 1; }
'''

JS = '''
export default function(component) {
    const { data, parentElement, setTriggerValue } = component;
    let root = parentElement.querySelector('.hp-root');
    if (!root) {
        root = document.createElement('div');
        root.className = 'hp-root';
        parentElement.appendChild(root);
    }
    if (root.dataset.markup !== data.markup) {
        root.innerHTML = data.markup;
        root.dataset.markup = data.markup;
    }
    const hp = root.querySelector('.hp');
    const chips = [...root.querySelectorAll('.hp-host')];
    const items = [...root.querySelectorAll('li[data-parasite]')];
    const hostsOf = (li) => (li.dataset.hosts || '').split('SEP').filter(Boolean);

    // the selection Python holds, redrawn on every run
    chips.forEach((c) => c.classList.toggle('sel', c.dataset.host === data.host));
    items.forEach((li) => {
        li.classList.toggle('sel', li.dataset.parasite === data.parasite);
        li.classList.toggle('off', !!data.host && !hostsOf(li).includes(data.host));
    });

    const clear = () => {
        hp.classList.remove('hovering-host', 'hovering-parasite');
        root.querySelectorAll('.lit').forEach((e) => e.classList.remove('lit'));
    };
    chips.forEach((chip) => {
        chip.onmouseenter = () => {
            hp.classList.add('hovering-host');
            chip.classList.add('lit');
            items.forEach((li) => li.classList.toggle('lit', hostsOf(li).includes(chip.dataset.host)));
        };
        chip.onmouseleave = clear;
        chip.onclick = () => setTriggerValue('pick', { host: chip.dataset.host });
    });
    items.forEach((li) => {
        li.onmouseenter = () => {
            // with a host picked the hosts stay as they are: the parasites are being
            // chosen among that host's
            if (li.classList.contains('na') || data.host) return;
            const hosts = hostsOf(li);
            hp.classList.add('hovering-parasite');
            chips.forEach((c) => c.classList.toggle('lit', hosts.includes(c.dataset.host)));
        };
        li.onmouseleave = clear;
        li.onclick = () => {
            if (li.classList.contains('na') || li.classList.contains('off')) return;
            setTriggerValue('pick', { parasite: li.dataset.parasite });
        };
    });
}
'''.replace('SEP', SEP)

def split_label(label):
    '''
    The common and the latin name of a host, from the label the configuration gives it.
    '''
    found = LABEL.match(label)
    if not found:
        return label, label

    return found['common'].capitalize(), found['latin']


def attr(value):
    return html.escape(value, quote=True)


def dots(host_ids, hosts, infected):
    '''
    A dot per host the species infects, in the order of the configuration.
    '''
    return ''.join(
        f'<span class="hp-dot" style="--c: {hosts[h]["color"]}" '
        f'title="{attr(split_label(hosts[h]["label"])[0])}"></span>'
        for h in host_ids if h in infected)


def host_card(host, count):
    common, latin = split_label(host['label'])
    return (f'<div class="hp-host" data-host="{attr(host["label"])}" style="--c: {host["color"]}">'
            f'<span class="hp-dot"></span><b>{html.escape(common)}</b>'
            f'<i>{html.escape(latin)}</i><span title="parasites">{count}</span></div>')


def parasite_item(name, infected, host_ids, hosts):
    '''
    One parasite of a group box; one without predictions is listed but cannot be picked.
    '''
    labels = SEP.join(hosts[h]['label'] for h in host_ids if h in infected)
    na = '' if infected else ' class="na" title="No predicted interactions"'
    return (f'<li data-parasite="{attr(name)}" data-hosts="{attr(labels)}"{na}>'
            f'{dots(host_ids, hosts, infected)}<i>{html.escape(name)}</i></li>')


def group_panel(group, color, members, host_ids, hosts):
    icon = f'<svg viewBox="0 0 100 100">{ICONS[GROUP_ICON[group]](color)}</svg>'
    items = ''.join(parasite_item(name, infected, host_ids, hosts) for name, infected in members)
    wide = ' wide' if group == WIDE_GROUP else ''
    return (f'<div class="hp-group{wide}" style="--c: {color}"><div class="hp-head">{icon}'
            f'<div><b>{html.escape(group)}</b><small>{html.escape(GROUP_SUB.get(group, ""))}</small></div>'
            f'<em>{len(members)}</em></div><ul>{items}</ul></div>')


@st.cache_data(show_spinner=False)
def predicted_hosts(data_dir):
    '''The host taxids each parasite has predictions against, by parasite label.'''
    predictions = web_utils.load_predictions(data_dir)
    pairs = predictions[['taxid1_label', 'taxid2']].drop_duplicates()

    return {parasite: {int(t) for t in rows['taxid2']}
            for parasite, rows in pairs.groupby('taxid1_label')}


def hosts_of(config, data_dir, parasite):
    '''The host labels a parasite has predictions against, in config order.'''
    infected = predicted_hosts(data_dir).get(parasite, set())

    return [host['label'] for taxid, host in config['hosts'].items() if taxid in infected]


def markup(config, data_dir):
    hosts = config['hosts']
    parasites = config['parasites']
    palette = config.get('parasite_groups', {})
    host_ids = list(hosts)
    infected_by = predicted_hosts(data_dir)

    by_group = {}
    for parasite in parasites.values():
        by_group.setdefault(parasite['group'], []).append(
            (parasite['label'], infected_by.get(parasite['label'], set())))
    # the groups the slide orders, then any the configuration has since gained
    groups = [g for g in GROUP_ORDER if g in by_group] + sorted(set(by_group) - set(GROUP_ORDER))

    cards = ''.join(
        host_card(host, sum(1 for infected in infected_by.values() if taxid in infected))
        for taxid, host in hosts.items())
    panels = {
        group: group_panel(group, palette.get(group, '#999999'), sorted(by_group[group]),
                           host_ids, hosts)
        for group in groups}
    rest = ''.join(panel for group, panel in panels.items() if group != WIDE_GROUP)

    return (f'<div class="hp"><div class="hp-label">HOSTS</div><div class="hp-hosts">{cards}</div>'
            f'<div class="hp-label">PARASITES</div><div class="hp-groups">'
            f'{panels.get(WIDE_GROUP, "")}<div class="hp-rest">{rest}</div></div></div>')


def read_pick():
    '''The click the component last sent, if any.'''
    result = st.session_state.get(PICKER_KEY)
    if result is None:
        return None
    if isinstance(result, dict):
        return result.get('pick')

    return getattr(result, 'pick', None)


def apply_pick(config, data_dir):
    '''
    Turns a click into the shared selection, before the page is drawn. A host click
    toggles the host and drops a parasite it does not have; a parasite click toggles the
    parasite and, when no host is chosen and the parasite has only one, chooses that host.
    '''
    pick = read_pick()
    if not pick:
        return

    host, parasite = web_utils.get_selection()
    if 'host' in pick:
        host = None if pick['host'] == host else pick['host']
        if host and parasite and host not in hosts_of(config, data_dir, parasite):
            parasite = None
    elif 'parasite' in pick:
        parasite = None if pick['parasite'] == parasite else pick['parasite']
        if parasite and host is None:
            infected = hosts_of(config, data_dir, parasite)
            if len(infected) == 1:
                host = infected[0]
    web_utils.set_selection(host, parasite)


def show(config, data_dir):
    '''Draws the overview, marking the current selection; clicks update it.'''
    host, parasite = web_utils.get_selection()
    # registered on every run, since a registration does not outlive the runtime that made
    # it; the same definition again is accepted silently
    component = st.components.v2.component('organism_picker', css=CSS, js=JS)
    component(key=PICKER_KEY,
              data={'markup': markup(config, data_dir), 'host': host, 'parasite': parasite},
              on_pick_change=lambda: apply_pick(config, data_dir))
