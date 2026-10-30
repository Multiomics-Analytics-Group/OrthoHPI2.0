'''
Draws a similarity matrix of parasites as a clickable heatmap, with strips of annotations
down the side and along the top, so that the two axes are visibly the same list in the same
order.
'''
import numpy as np
import plotly.graph_objects as go
import streamlit as st

# side of a heatmap cell in pixels, which sizes the figure
CELL = 22
# the diagonal of the heatmap: a grey, so it is not read as a count
DIAGONAL_COLOUR = '#9e9e9e'
# pixels the colour bar takes right of the matrix; held in pixels so plotly does not widen
# the margin and unsquare the plot
COLORBAR_GAP = 8
COLORBAR_ROOM = 71
COLORBAR_DIGIT = 6.5
# two decimals: most pairs sit under a tenth
TICK = '.2f'
# height of one legend row below the heatmap, in pixels
LEGEND_ROW = 22
# room one legend entry takes, in pixels: marker, padding and ~7.5 px a character of the
# longest name
LEGEND_ENTRY = 40
LEGEND_CHAR = 7.5
# smallest the names on the column axis are written, in points
SMALLEST_LABEL = 9


@st.cache_data(show_spinner=False)
def flat_colour_scale(values, colors):
    '''
    The colour scale of a strip of flat bands -- one band per value, nothing
    interpolated between two of them -- and the codes that index it, which is what a
    strip is drawn from.
    '''
    steps = [step for i, value in enumerate(colors)
             for step in ([i / len(colors), colors[value]],
                          [(i + 1) / len(colors), colors[value]])]

    return dict(colorscale=steps, zmin=-0.5, zmax=len(colors) - 0.5, showscale=False,
                hovertemplate='%{text}<extra></extra>'), \
        [list(colors).index(v) for v in values]


def short_parasite(label):
    '''`Plasmodium falciparum` as `P. falciparum`, the way a column is named.'''
    return f'{label[0]}. {label.split(" ")[1]}'


def generate(similarity, shared, strips, column, x_names=None,
             count_label='Shared interactors', max_side=None, smallest_label=SMALLEST_LABEL,
             x_angle=-60):
    '''
    The similarity matrix, with a strip per annotation in `strips` -- the value of each
    row and the colours of the values, in legend order -- beside the rows and above the
    columns, the matrix no wider than `column` nor than `max_side` pixels and the column
    names slanted `x_angle` degrees. Returns the figure and the index of the layer a click
    on a cell comes from.
    '''
    strips = [(values, {v: c for v, c in colors.items() if v in set(values)})
              for values, colors in strips]
    y_names = list(similarity.index)
    x_names = x_names or [short_parasite(g) for g in y_names]
    cells = list(range(len(y_names)))

    # the same span on both axes, so a square plot area is one of square cells
    span = len(cells) + 0.6 + len(strips)
    # the ticks beside the colour bar are all that varies in the room it needs
    right = COLORBAR_ROOM + COLORBAR_DIGIT * len(f'{np.nanmax(similarity.to_numpy()):{TICK}}')

    # a line of text is about 1.1 times its point size tall
    def fits(room, smallest=6):
        sizes = [size for size in (10, 9, 8, 7, 6) if size >= smallest]
        return next(size for size in sizes if room >= 1.1 * size or size == sizes[-1])

    longest = max(len(name) for name in y_names)

    def square(size):
        '''Side of the square and the left margin, with the names written `size` points.'''
        # what the cells are held square to
        margin = 0.65 * size * longest + 12

        return max(240, min(CELL * span, column - margin - right,
                            max_side or float('inf'))), margin

    # measured once at the largest size and again at the size that left them
    font = fits(square(10)[0] / span)
    side, left = square(font)
    # column names are held to a floor of their own, being the only thing that names a
    # column
    x_font = fits(side / span, smallest=smallest_label)

    # slanted, the names stand the sine of their angle of their length tall
    top = abs(np.sin(np.radians(x_angle))) * (
        0.65 * x_font * max(len(name) for name in x_names) + 12) + 25
    # plotly wraps legend entries to the plot width, so the rows are counted rather than
    # assumed
    def legend_rows(names):
        entry = LEGEND_ENTRY + LEGEND_CHAR * max(len(name) for name in names)

        return -(-len(names) // max(1, int(side // entry)))

    rows = [legend_rows(list(colors)) for _, colors in strips]
    bottom = 34 + LEGEND_ROW * sum(rows)

    figure = go.Figure()
    # a strip beside the rows and above the columns per annotation, a cell clear of the
    # matrix
    for offset, (values, colors) in enumerate(strips):
        scale, codes = flat_colour_scale(values, colors)
        figure.add_trace(go.Heatmap(z=[[code] for code in codes],
                                    x0=-1.4 - offset, dx=1, y=cells,
                                    text=[[v] for v in values], ygap=1, **scale))
        figure.add_trace(go.Heatmap(z=[codes], x=cells,
                                    y0=-1.4 - offset, dy=1,
                                    text=[list(values)], xgap=1, **scale))

    # hover belongs to the clickable layer below. The scale runs to the largest similarity
    # rather than 1.0
    figure.add_trace(go.Heatmap(z=similarity.to_numpy(), x=cells, y=cells, hoverinfo='skip',
                                colorscale=['#ffffff', '#deebf7', '#9ecae1', '#6baed6',
                                            '#3182bd', '#08519c'],
                                zmin=0, zmax=np.nanmax(similarity.to_numpy()),
                                hoverongaps=False,
                                colorbar=dict(title=dict(text='Jaccard similarity',
                                                         side='right'),
                                              thickness=12, len=0.6, y=1, yanchor='top',
                                              x=1 + COLORBAR_GAP / side,
                                              tickformat=TICK,
                                              tickfont=dict(size=10))))

    # the diagonal, drawn over the empty cells; both traces turn hovering on gaps off so
    # neither answers for the cells beneath
    diagonal = np.full(similarity.shape, np.nan)
    np.fill_diagonal(diagonal, 0)
    figure.add_trace(go.Heatmap(z=diagonal, x=cells, y=cells,
                                text=[[name] * len(y_names) for name in y_names],
                                colorscale=[[0, DIAGONAL_COLOUR], [1, DIAGONAL_COLOUR]],
                                zmin=0, zmax=1, showscale=False, hoverongaps=False,
                                hovertemplate='%{text}<extra></extra>'))

    # a heatmap cannot be clicked in Streamlit, so every cell carries an invisible scatter
    # marker; x and y are the indices of the two parasites
    click_targets = [(x, y) for y in cells for x in cells if x != y]
    figure.add_trace(go.Scatter(
        x=[x for x, _ in click_targets], y=[y for _, y in click_targets], mode='markers',
        marker=dict(symbol='square', size=side / span, color='rgba(0,0,0,0)',
                    line=dict(width=0)),
        # plotly would otherwise dim everything but the clicked cell
        selected=dict(marker=dict(opacity=1)), unselected=dict(marker=dict(opacity=1)),
        # the count beside the ratio, which is what the dialog then lists
        text=[f'{y_names[y]} and {y_names[x]}<br>'
              f'Jaccard similarity {similarity.iat[y, x]:{TICK}}<br>'
              f'{count_label} {shared.iat[y, x]:.0f}' for x, y in click_targets],
        hovertemplate='%{text}<extra></extra>', showlegend=False))
    click_layer = len(figure.data) - 1

    # the strips are heatmaps and cannot carry a legend, so their values are named by empty
    # traces, one legend per strip
    legends = ['legend'] + [f'legend{i}' for i in range(2, len(strips) + 1)]
    for legend, (_, colors) in zip(legends, strips):
        for value, color in colors.items():
            figure.add_trace(go.Scatter(x=[None], y=[None], mode='markers', name=value,
                                        marker=dict(size=10, symbol='square', color=color),
                                        legend=legend, hoverinfo='skip', showlegend=True))

    # every parasite is named, however small the cells
    ticks = dict(tickmode='array', tickvals=cells, ticks='')
    start = -1 - len(strips)
    # named above the matrix, slanted the way the other figures of the page slant them
    figure.update_xaxes(range=[start, len(cells) - 0.4], side='top',
                        ticktext=x_names, tickangle=x_angle, tickfont=dict(size=x_font), **ticks)
    # reversed, so the first parasite is the top row; the range runs back past the strips
    figure.update_yaxes(range=[len(cells) - 0.4, start],
                        ticktext=y_names, tickfont=dict(size=font), **ticks)

    # the plot area is `side` pixels each way, the margins holding everything else, which
    # keeps the cells square; each legend sits under the rows of the ones before it
    entries = dict(orientation='h', yanchor='top', xanchor='left', x=0, itemclick=False,
                   itemdoubleclick=False, font=dict(size=11), title_font=dict(size=11))
    placed = {legend: dict(y=-0.01 - LEGEND_ROW * sum(rows[:i]) / side, **entries)
              for i, legend in enumerate(legends)}
    figure.update_layout(width=left + side + right, height=side + top + bottom,
                         plot_bgcolor='white',
                         # no pad, so the names sit against the strips
                         margin=dict(l=left, r=right, t=top, b=bottom, pad=0),
                         **placed)

    return figure, click_layer
