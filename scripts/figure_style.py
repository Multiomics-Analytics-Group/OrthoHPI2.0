'''
What the paper's figures share: the sizes Frontiers asks for (text no smaller than 8 pt
at final size, 180 mm across two columns, nothing longer than a page), the row height
that gives 8 pt names room, and the abbreviated species name the figures read under.
'''
import matplotlib

# Frontiers: two columns are 180 mm, and a figure is at most a page
WIDTH = 180 / 25.4
MAX_HEIGHT = 225 / 25.4
# the smallest text Frontiers allows, and one step up for titles
FONT = 8
TITLE = 9
# inches per row of 8 pt names, with a gap between them
ROW = 0.17


def apply():
    matplotlib.rcParams.update({
        'font.family': 'sans-serif',
        'font.size': FONT,
        'axes.titlesize': TITLE,
        'axes.labelsize': FONT,
        'xtick.labelsize': FONT,
        'ytick.labelsize': FONT,
        'legend.fontsize': FONT,
        'hatch.linewidth': 0.4,
    })


def short_name(species):
    '''`Trichinella spiralis` as `T. spiralis`, the abbreviation the app uses.'''
    parts = species.split(' ')

    return f'{parts[0][0]}. {" ".join(parts[1:])}' if len(parts) > 1 else species


def stack_legends(fig, blocks, ncol=4, left=0.1, bottom=0.02):
    '''
    One legend per kind of colour, one under another at the foot of the figure and
    flush with the same left edge (`left`, a figure fraction), so that no two kinds
    share a column. `blocks` is [(title, handles)] or [(title, handles, ncol)] for a
    block with entries too long for the default columns, top block first. Returns the
    height taken, in inches, for the caller to leave under its axes.
    '''
    height = fig.get_figheight()
    line = 1.1 * FONT / 72
    y = bottom
    for block in reversed(blocks):
        title, handles = block[:2]
        columns = min(block[2] if len(block) > 2 else ncol, len(handles))
        rows = -(-len(handles) // columns)
        fig.legend(handles=handles, title=title, loc='lower left', ncol=columns,
                   frameon=False, handlelength=0.9, columnspacing=1.2, alignment='left',
                   bbox_to_anchor=(left, y / height), title_fontsize=FONT)
        y += (rows + 1) * line + 0.12

    return y
