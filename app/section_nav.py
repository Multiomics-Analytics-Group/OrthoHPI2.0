'''
The layout of a long analysis page: a narrow column on the left that stays in view while
the page scrolls, listing the sections of the page and marking the one being read, beside
the column the figures are drawn in.
'''
import streamlit as st

import web_utils

# the share of the page the two columns take
NAV_SHARE, BODY_SHARE = 1, 4.3
# what streamlit leaves between columns at gap='medium', in pixels
NAV_GAP = 32
# a section counts as being read once its heading is this far from the top of the window,
# below the bar across the top
READ_LINE = 140

# drawn above every section heading but the first
DIVIDER = ('<div style="border-top: 1px solid rgba(128, 128, 128, 0.2); '
           'margin: 1.5rem 0 0.25rem;"></div>')

CSS = '''
nav { display: flex; flex-direction: column; gap: 0.1rem; font-size: 0.9rem; }
.label { font-size: 0.75rem; font-weight: 600; letter-spacing: 0.1em; opacity: 0.6;
         margin-bottom: 0.3rem; }
a { color: inherit; text-decoration: none; padding: 0.3rem 0.6rem; border-radius: 6px;
    border-left: 3px solid transparent; line-height: 1.3; cursor: pointer; opacity: 0.75; }
a:hover { background: color-mix(in srgb, currentColor 6%, transparent); opacity: 1; }
a.active { border-left-color: var(--st-primary-color); opacity: 1; font-weight: 600;
           background: color-mix(in srgb, var(--st-primary-color) 8%, transparent); }
'''

JS = '''
export default function(component) {
    const { data, parentElement } = component;
    let nav = parentElement.querySelector('nav');
    if (!nav) {
        nav = document.createElement('nav');
        parentElement.appendChild(nav);
    }
    nav.innerHTML = '<div class="label">ON THIS PAGE</div>' + data.sections.map(
        ([title, anchor]) => `<a data-target="${anchor}">${title}</a>`).join('');
    const links = [...nav.querySelectorAll('a')];

    links.forEach((link) => {
        link.onclick = () => {
            const target = document.getElementById(link.dataset.target);
            if (!target) return;
            target.style.scrollMarginTop = '5rem';
            // charts above can still be resizing, which moves the heading after a jump,
            // so it is jumped to again once they have settled
            target.scrollIntoView({ block: 'start' });
            setTimeout(() => target.scrollIntoView({ block: 'start' }), 400);
        };
    });

    // the section being read is the last one whose heading has passed the read line
    const mark = () => {
        let current = links[0];
        for (const link of links) {
            const target = document.getElementById(link.dataset.target);
            if (target && target.getBoundingClientRect().top < READ_LINE) current = link;
        }
        links.forEach((link) => link.classList.toggle('active', link === current));
    };
    mark();
    // scroll events do not bubble, so they are caught on the way down to the scrolling
    // element, whichever it is
    document.addEventListener('scroll', mark, true);

    return () => document.removeEventListener('scroll', mark, true);
}
'''.replace('READ_LINE', str(READ_LINE))


class SectionNav:
    '''
    Splits the page into the section list and the body. The list is drawn first, from the
    sections the page drew last time, or `expected` on a first visit, so it is in place
    before the figures; sections are named through `subheader` as the body is drawn, and
    `finish` redraws the list only where they turned out different.
    '''

    def __init__(self, page, expected):
        nav, self.body = st.columns([NAV_SHARE, BODY_SHARE], gap='medium')
        column = nav.container(key='section_nav')
        # the way back to the start page, above the list
        self.back = column.container(key='section_back')
        self.links = column.empty()
        # the controls that set the page, at the foot of the column, apart from the list
        self.side = column.container(key='section_controls')
        self.state_key = f'sections_{page}'
        self.shown = [tuple(s) for s in st.session_state.get(self.state_key, expected)]
        self.draw(self.shown)
        self.sections = []

    def draw(self, sections):
        component = st.components.v2.component('section_nav', css=CSS, js=JS)
        # keyed by the list, so the same list is left mounted from one run to the next and
        # a changed one can replace it within a run
        key = 'section_nav_' + '_'.join(anchor for _, anchor in sections)
        with self.links:
            component(key=key, data={'sections': sections})

    @staticmethod
    def content_width():
        '''Width in pixels of the body column, which the figures are sized to.'''
        share = BODY_SHARE / (NAV_SHARE + BODY_SHARE)

        return (web_utils.page_width() - NAV_GAP) * share

    def subheader(self, title, anchor, help=None, slot=None, label=None):
        '''
        A section heading, listed in the column on the left under `label`, or under the
        title where that carries no markdown.
        '''
        # a placeholder holds one element, so the line and the heading share a box in it
        box = slot.container() if slot is not None else st
        if self.sections:
            self.divider(box)
        box.subheader(title, anchor=anchor, help=help)
        self.add(label or title, anchor)

    @staticmethod
    def divider(box=st):
        '''The pale line that parts a section from the one before it.'''
        box.markdown(DIVIDER, unsafe_allow_html=True)

    def add(self, label, anchor):
        '''Lists a heading drawn elsewhere, under the anchor it was given.'''
        if anchor not in (a for _, a in self.sections):
            self.sections.append((label, anchor))

    def finish(self):
        if self.sections != self.shown:
            st.session_state[self.state_key] = self.sections
            self.draw(self.sections)
