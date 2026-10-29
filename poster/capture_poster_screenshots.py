'''
Screenshots of the network page for the poster's web resource panels, written to
poster/figures/web-pf/: P. falciparum in human at confidence 0.7, restricted to blood and
to erythrocytes, as in the paper's case study. The AlphaFold dialog is opened for the
SRA - CUL4A link, a merozoite surface antigen in a small
component of the network, so the tooltip and the models are of one clean link. Reuses scripts/capture_web_screenshots.py, which takes
the paper's T. cruzi screenshots.

  ~/.venvs/streamlit/bin/streamlit run app/OrthoHPI_Home.py --server.headless true --server.port 8599 --theme.base light
  ~/.venvs/streamlit/bin/python poster/capture_poster_screenshots.py
'''
import os
import sys

from PIL import Image
from playwright.sync_api import TimeoutError as PlaywrightTimeout, sync_playwright

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.capture_web_screenshots import (DIALOG_WIDTH, FIND_LINK, NETWORK_PAGE, VIEWPORT,
                                             choose, collapse, expose_network,
                                             hide_captions, slide, wait)

URL = 'http://localhost:8599'
# a window this wide gives the page's column about the shape of the paper's screenshots
WINDOW = {**VIEWPORT, 'width': 860}
OUT = os.path.join('poster', 'figures', 'web-pf')
# PF3D7_1431400 (Surface-related antigen SRA) - CUL4A, confidence 0.943
LINK = ['36329.Q8ILF1', '9606.ENSP00000364589']


def clip(page, top, bottom, pad=12):
    '''The page's content column, as wide as the network, from one element to another.'''
    column = page.locator('iframe[src*="ppi_network"]').bounding_box()
    t, b = top.bounding_box(), bottom.bounding_box()
    return {'x': column['x'] - pad, 'y': t['y'] - pad, 'width': column['width'] + 2 * pad,
            'height': b['y'] + b['height'] + 2 * pad - t['y']}


def pick(page, index, text):
    '''Types `text` into the `index`th multiselect of the page and takes the match.'''
    page.locator('[data-testid=stMultiSelect]').nth(index).locator('input').click()
    page.keyboard.type(text)
    page.wait_for_timeout(500)
    page.keyboard.press('Enter')
    page.keyboard.press('Escape')
    wait(page)


def network_page(page):
    page.goto(f'{URL}/Predicted_Host-Parasite_PPIs')
    wait(page)
    choose(page, 'Select a host', 'Human')
    choose(page, 'Select a parasite', 'falciparum')
    slide(page, 0, 35)   # confidence 0.35 -> 0.70
    pick(page, 0, 'blood')
    pick(page, 1, 'erythro')

    frame = next(f for f in page.frames if 'ppi_network' in f.url)
    network = page.locator('iframe[src*="ppi_network"]')
    link_ids = LINK

    dialog = page.get_by_role('dialog')
    for _ in range(5):
        page.wait_for_timeout(4000)  # the network fits itself once the physics settle
        link = frame.evaluate(FIND_LINK, link_ids)
        box = network.bounding_box()
        page.mouse.click(box['x'] + link['x'], box['y'] + link['y'])
        try:
            dialog.wait_for(timeout=30000)
            break
        except PlaywrightTimeout:
            pass
    wait(page)
    page.set_viewport_size({**VIEWPORT, 'width': DIALOG_WIDTH})
    hide_captions(page)
    page.evaluate('document.activeElement.blur()')
    page.wait_for_timeout(5000)  # 3Dmol draws after the run ends, and again on the resize
    names = dialog.locator('h4').first.bounding_box()
    buttons = dialog.get_by_text('PDB file').first.bounding_box()
    frame_box = dialog.bounding_box()
    page.screenshot(path=f'{OUT}/alphafold.png', clip={
        'x': frame_box['x'], 'y': names['y'] - 6, 'width': frame_box['width'],
        'height': buttons['y'] - names['y'] - 6})
    collapse(Image.open(f'{OUT}/alphafold.png')).save(f'{OUT}/alphafold.png')
    page.keyboard.press('Escape')
    page.set_viewport_size(WINDOW)
    wait(page)

    # hover the link so its tooltip shows in the network screenshot
    tooltip = frame.locator('.vis-tooltip')
    for _ in range(5):
        page.wait_for_timeout(4000)
        box = network.bounding_box()
        link = frame.evaluate(FIND_LINK, link_ids)
        page.mouse.move(0, 0)
        page.mouse.move(box['x'] + link['x'], box['y'] + link['y'], steps=5)
        page.wait_for_timeout(1500)
        if tooltip.count() and tooltip.first.is_visible():
            break

    # down to the network's heading, so the body figure's colour scale is in
    heading = page.get_by_role('heading', name='Network of host-parasite PPIs')
    controls = clip(page, page.get_by_text('Predicted interactions per infected organ'), heading)
    controls['height'] -= heading.bounding_box()['height'] + 12
    page.screenshot(path=f'{OUT}/controls.png', clip=controls)
    collapse(Image.open(f'{OUT}/controls.png'), margin=40).save(f'{OUT}/controls.png')
    page.screenshot(path=f'{OUT}/network.png', clip=clip(
        page, page.get_by_text('Extracellular', exact=True).first, network))
    collapse(Image.open(f'{OUT}/network.png'), margin=40).save(f'{OUT}/network.png')


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page(viewport=WINDOW, device_scale_factor=2)
        page.route(NETWORK_PAGE, expose_network)
        network_page(page)
        browser.close()
    print(f'Wrote the screenshots to {OUT}')
