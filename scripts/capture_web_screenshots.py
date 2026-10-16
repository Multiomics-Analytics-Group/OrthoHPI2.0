'''
Screenshots of the network page for the web interface figure, written to
paper/figures/web/ for scripts/build_web_interface_figure.py: T. cruzi in human at
confidence 0.7, restricted to the heart.

  controls    the body figure and the query controls
  network     the network, the mouse over the clicked link for its tooltip
  alphafold   the AlphaFold models of the Q4D3L1 - DR1 link, opened by clicking it: the
              name of each protein over its model
  enrichment  the chart of the GO biological processes enriched among the host proteins

Needs the app running and Playwright with Chromium, both in the streamlit environment:

  streamlit run app/OrthoHPI_Home.py --server.headless true --server.port 8599 --theme.base light
  ~/.venvs/streamlit/bin/python scripts/capture_web_screenshots.py
'''
import argparse
import os
import re

import numpy as np
from PIL import Image
from playwright.sync_api import TimeoutError as PlaywrightTimeout, sync_playwright

VIEWPORT = {'width': 1500, 'height': 9000}
# the window width the AlphaFold dialog is captured at, narrowing it to half the figure's width
DIALOG_WIDTH = 900
# the centred column of the network page, in CSS pixels
COLUMN = (390, 1110)
HIDE_CAPTIONS = '[data-testid=stCaptionContainer]{display:none !important}'
# the network component keeps its vis.js network inside a closure; the served page is
# patched to hand it out, so a link can be found and clicked
NETWORK_PAGE = re.compile(r'/component/ppi_network[^/]*/index\.html')
# the link clicked for its AlphaFold models: the best-modelled link of the structural benchmark
LINK = ('353153.Q4D3L1', '9606.ENSP00000359295')
# a point on the link between two nodes that the network takes for that link and no node,
# as a click there is read; failing such a link, the point of the link whose midpoint lies
# furthest from any node
FIND_LINK = '''([u, v]) => {
  const net = window.__getNet();
  net.stopSimulation();
  const pos = net.getPositions();
  const box = document.getElementById('mynetwork').getBoundingClientRect();
  const dom = id => net.canvasToDOM(pos[id]);
  const edge = net.body.data.edges.get().find(e => (e.from === u && e.to === v) || (e.from === v && e.to === u));
  if (edge) {
    const a = dom(u), b = dom(v);
    for (const t of [0.5, 0.4, 0.6, 0.3, 0.7, 0.25, 0.75]) {
      const p = {x: a.x + t * (b.x - a.x), y: a.y + t * (b.y - a.y)};
      if (net.getNodeAt(p) === undefined && net.getEdgeAt(p) === edge.id) {
        return {x: p.x + box.left, y: p.y + box.top};
      }
    }
  }
  let best = null;
  for (const e of net.body.data.edges.get()) {
    const a = dom(e.from), b = dom(e.to);
    const m = {x: (a.x + b.x) / 2, y: (a.y + b.y) / 2};
    let near = Math.hypot(a.x - b.x, a.y - b.y) / 2;
    for (const id in pos) { const p = dom(id); near = Math.min(near, Math.hypot(p.x - m.x, p.y - m.y)); }
    if (!best || near > best.near) best = {near, x: m.x + box.left, y: m.y + box.top};
  }
  return best;
}'''


def expose_network(route):
    response = route.fetch()
    body = response.text().replace(
        'var network = null;', 'var network = null; window.__getNet = function () { return network; };', 1)
    route.fulfill(response=response, body=body)


def wait(page):
    '''Until Streamlit has finished the run.'''
    page.wait_for_timeout(1000)
    page.wait_for_function("() => !document.querySelector('[data-testid=stStatusWidget]')", timeout=240000)
    page.wait_for_timeout(1000)


def choose(page, label, text):
    page.get_by_label(label).first.click()
    page.keyboard.type(text)
    page.wait_for_timeout(300)
    page.keyboard.press('Enter')
    wait(page)


def slide(page, index, steps):
    page.get_by_role('slider').nth(index).focus()
    for _ in range(steps):
        page.keyboard.press('ArrowRight')
    wait(page)


def hide_captions(page):
    page.add_style_tag(content=HIDE_CAPTIONS)
    page.wait_for_timeout(1500)


def column_clip(top, bottom, pad=12):
    '''The network page's column from the top of one element to the bottom of another.'''
    t, b = top.bounding_box(), bottom.bounding_box()
    return {'x': COLUMN[0], 'y': t['y'] - pad, 'width': COLUMN[1] - COLUMN[0],
            'height': b['y'] + b['height'] + 2 * pad - t['y']}


def collapse(image, keep=60, margin=60, tol=6):
    '''Shrinks each band of blank rows taller than `keep` to `keep`, testing all but the
    `margin` at either side.'''
    pixels = np.asarray(image.convert('RGB')).astype(int)
    core = pixels[:, margin:pixels.shape[1] - margin]
    blank = (core.max(axis=1) - core.min(axis=1)).max(axis=1) <= tol
    rows, run = [], []
    for i, is_blank in enumerate(blank):
        if is_blank:
            run.append(i)
            continue
        rows += run[:keep // 2] + run[-(keep // 2):] if len(run) > keep else run
        rows.append(i)
        run = []
    rows += run[:keep]

    return Image.fromarray(pixels[rows].astype(np.uint8))


def network_page(page, base, out):
    page.goto(f'{base}/Predicted_Host-Parasite_PPIs')
    wait(page)
    choose(page, 'Select a host', 'Human')
    choose(page, 'Select a parasite', 'cruzi')
    slide(page, 0, 35)   # confidence 0.35 -> 0.70
    page.get_by_label('Select tissues to filter').first.click()
    page.keyboard.type('heart')
    page.wait_for_timeout(300)
    page.keyboard.press('Enter')
    page.keyboard.press('Escape')
    wait(page)

    # the AlphaFold dialog
    frame = next(f for f in page.frames if 'ppi_network' in f.url)
    network = page.locator('iframe[src*="ppi_network"]')
    dialog = page.get_by_role('dialog')
    for _ in range(5):
        page.wait_for_timeout(4000)  # the network fits itself once the physics settle
        link = frame.evaluate(FIND_LINK, list(LINK))
        box = network.bounding_box()
        page.mouse.click(box['x'] + link['x'], box['y'] + link['y'])
        try:
            dialog.wait_for(timeout=30000)
            break
        except PlaywrightTimeout:
            pass
    dialog.wait_for(timeout=1000)
    wait(page)
    page.set_viewport_size({**VIEWPORT, 'width': DIALOG_WIDTH})
    hide_captions(page)
    page.evaluate('document.activeElement.blur()')
    page.wait_for_timeout(5000)  # 3Dmol draws after the run ends, and again on the resize
    # from the protein names down to the buttons under the models
    names = dialog.locator('h4').first.bounding_box()
    buttons = dialog.get_by_text('PDB file').first.bounding_box()
    frame_box = dialog.bounding_box()
    page.screenshot(path=f'{out}/alphafold.png', clip={
        'x': frame_box['x'], 'y': names['y'] - 6, 'width': frame_box['width'],
        'height': buttons['y'] - names['y'] - 6})
    collapse(Image.open(f'{out}/alphafold.png')).save(f'{out}/alphafold.png')
    page.keyboard.press('Escape')
    page.set_viewport_size(VIEWPORT)
    wait(page)
    # off the link and back onto it, so the network sees a hover and shows the link's tooltip
    tooltip = frame.locator('.vis-tooltip')
    for _ in range(5):
        page.wait_for_timeout(4000)  # the network fits itself again after the resize
        box = network.bounding_box()
        link = frame.evaluate(FIND_LINK, list(LINK))
        page.mouse.move(0, 0)
        page.mouse.move(box['x'] + link['x'], box['y'] + link['y'], steps=5)
        page.wait_for_timeout(1500)
        if tooltip.count() and tooltip.first.is_visible():
            break

    page.screenshot(path=f'{out}/controls.png', clip=column_clip(
        page.get_by_text('Predicted interactions per infected organ'),
        page.locator('[data-testid=stMultiSelect]', has=page.get_by_label('Select cell types to filter'))))
    page.screenshot(path=f'{out}/network.png', clip=column_clip(
        page.get_by_text('Extracellular', exact=True).first, network))
    collapse(Image.open(f'{out}/network.png'), margin=40).save(f'{out}/network.png')
    heading = page.get_by_role('heading', name='Enriched biological processes', exact=True).bounding_box()
    # the first chart under the heading, the ranked processes of the host proteins, alone:
    # the panel title and the caption say what the heading and tabs above it do
    chart = next(box for box in (c.bounding_box() for c in page.locator('[data-testid=stPlotlyChart]').all())
                 if box and box['y'] > heading['y'])
    top = chart['y'] - 12
    bottom = chart['y'] + chart['height'] + 12
    page.screenshot(path=f'{out}/enrichment.png', clip={'x': COLUMN[0], 'y': top,
                                               'width': COLUMN[1] - COLUMN[0], 'height': bottom - top})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--url', default='http://localhost:8599')
    parser.add_argument('--out', default='paper/figures/web')
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page(viewport=VIEWPORT, device_scale_factor=2)
        page.route(NETWORK_PAGE, expose_network)
        network_page(page, args.url, args.out)
        browser.close()
    print(f'Wrote the screenshots to {args.out}')
