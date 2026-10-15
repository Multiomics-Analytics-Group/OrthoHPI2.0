'''
The app's frontal human body figure with the liver shaded, for the tissue step of the
workflow figure. Writes paper/figures/workflow_body.svg, which
paper/figures/workflow.svg links to.
'''
import os
import sys
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'app'))

import body_figure

OUT = os.path.join(ROOT, 'paper', 'figures', 'workflow_body.svg')
# the host colour of the workflow figure
SHADE = '#0072B2'


def main():
    # load_figure reads its figures relative to the repo root, as the app does
    os.chdir(ROOT)
    svg, _ = body_figure.load_figure('human')
    root = ET.fromstring(body_figure.shade_figure(svg, {'liver': 1}, [(1, SHADE)]))

    _, y, _, height = (float(v) for v in root.get('viewBox').replace(',', ' ').split())
    x, width = body_figure.HUMAN_FRONT_VIEW
    root.set('viewBox', f'{x} {y} {width} {height}')
    root.attrib.pop('style', None)
    root.set('width', str(width))
    root.set('height', str(height))

    ET.ElementTree(root).write(OUT, encoding='utf-8', xml_declaration=True)
    print(f'wrote {OUT} ({width:.0f} x {height:.0f})')


if __name__ == '__main__':
    main()
