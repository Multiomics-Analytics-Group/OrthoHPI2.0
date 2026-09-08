"""
Build a slide figure of the four OrthoHPI 2.0 hosts and all their parasites,
with the parasites grouped by taxonomy.

Everything is read from config.yml, so the figure follows the study whenever
hosts or parasites are added.  The organism drawings are simple hand-drawn
vector icons (no external assets), so the output is a single self-contained
SVG that PowerPoint / Keynote / Illustrator can open and recolour.

    .venv/bin/python scripts/build_hosts_parasites_figure.py [-o docs/hosts_parasites.svg]
"""
import argparse
import functools
import math
import os
import re

import yaml

# --- style -------------------------------------------------------------------
W, H = 1600, 900
MARGIN = 56
INK = "#22262b"
MUTED = "#6b737c"
PANEL_BG = "#f6f7f9"
PANEL_LINE = "#e2e5e9"
FONT = "'Helvetica Neue', Helvetica, Arial, sans-serif"

HOST_ORDER = [9606, 10116, 10090, 9823]
COMMON_NAME = {9606: "Human", 10116: "Rat", 10090: "Mouse", 9823: "Pig"}
# the host colours a parasite row can be marked with, one dot per host
DOT_HOSTS = [(9606, "Human"), (10116, "Rat"), (10090, "Mouse"), (9823, "Pig")]

GROUP_ORDER = ["Nematoda", "Trematoda", "Cestoda", "Apicomplexa",
               "Kinetoplastida", "Other protozoa", "Microsporidia"]
# which icon to draw for each taxonomic group
GROUP_ICON = {
    "Nematoda": "nematode", "Trematoda": "fluke", "Cestoda": "tapeworm",
    "Apicomplexa": "apicomplexan", "Kinetoplastida": "trypanosome",
    "Other protozoa": "amoeba", "Microsporidia": "spore",
}
GROUP_SUB = {
    "Nematoda": "roundworms", "Trematoda": "flukes", "Cestoda": "tapeworms",
    "Apicomplexa": "apicomplexans", "Kinetoplastida": "trypanosomatids",
    "Other protozoa": "amoeba, flagellates", "Microsporidia": "spore-forming fungi",
}


# --- little vector helpers ---------------------------------------------------
def ribbon(centreline, width, n=64):
    """Closed path for a tapered body: centreline(t) -> (x, y), width(t) -> half-width."""
    left, right = [], []
    for i in range(n + 1):
        t = i / n
        x, y = centreline(t)
        dt = 1e-3
        x1, y1 = centreline(min(1.0, t + dt))
        x0, y0 = centreline(max(0.0, t - dt))
        dx, dy = x1 - x0, y1 - y0
        norm = math.hypot(dx, dy) or 1.0
        nx, ny = -dy / norm, dx / norm
        w = width(t)
        left.append((x + nx * w, y + ny * w))
        right.append((x - nx * w, y - ny * w))
    pts = left + right[::-1]
    d = "M{:.2f},{:.2f} ".format(*pts[0])
    d += " ".join("L{:.2f},{:.2f}".format(x, y) for x, y in pts[1:])
    return d + " Z"


def taper(peak=1.0, head=0.35, tail=0.08):
    """Half-width profile: thickest around the middle, blunt head, fine tail."""
    def f(t):
        base = math.sin(math.pi * (0.15 + 0.85 * t)) ** 0.7
        ends = head + (1 - head) * min(1.0, t / 0.18) if t < 0.18 else 1.0
        tail_f = tail + (1 - tail) * min(1.0, (1 - t) / 0.35)
        return peak * base * ends * tail_f
    return f


# --- organism icons ----------------------------------------------------------
# Every icon draws inside a 100 x 100 box and is placed with a transform.
def icon_nematode(c):
    body = ribbon(lambda t: (8 + 84 * t, 50 + 22 * math.sin(2 * math.pi * 1.1 * t + 0.5)),
                  taper(peak=10.5, head=0.75, tail=0.06))
    return (f'<path d="{body}" fill="{c}"/>'
            f'<circle cx="12" cy="59" r="2.6" fill="#ffffff" opacity="0.85"/>')


def icon_fluke(c):
    body = ("M50,6 C60,10 70,22 74,40 C80,64 72,88 50,95 "
            "C28,88 20,64 26,40 C30,22 40,10 50,6 Z")
    return (f'<path d="{body}" fill="{c}"/>'
            f'<circle cx="50" cy="19" r="5.5" fill="none" stroke="#ffffff"'
            f' stroke-width="3.2" opacity="0.95"/>'
            f'<circle cx="50" cy="44" r="8" fill="none" stroke="#ffffff"'
            f' stroke-width="3.4" opacity="0.95"/>'
            f'<path d="M50,56 C40,66 37,80 40,92 M50,56 C60,66 63,80 60,92"'
            f' stroke="#ffffff" stroke-width="3" fill="none" opacity="0.6"'
            f' stroke-linecap="round"/>')


def icon_tapeworm(c):
    parts = [f'<circle cx="13" cy="50" r="9" fill="{c}"/>',
             f'<circle cx="10" cy="45" r="2.6" fill="#ffffff" opacity="0.9"/>',
             f'<circle cx="10" cy="55" r="2.6" fill="#ffffff" opacity="0.9"/>',
             f'<path d="M20,47 h6 v6 h-6 z" fill="{c}"/>']
    x = 26
    for i in range(6):                      # proglottids widening down the strobila
        h = 12 + i * 4.6
        w = 9 + i * 1.6
        y = 50 - h / 2 + 3 * math.sin(i * 0.9)
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}"'
                     f' rx="3" fill="{c}"/>')
        x += w + 2.2
    return "".join(parts)


def icon_apicomplexan(c):
    body = ("M22,86 C10,64 14,32 38,14 C46,8 54,10 56,16 C58,22 52,26 46,32 "
            "C30,48 26,68 30,86 C31,92 24,92 22,86 Z")
    return (f'<path d="{body}" fill="{c}"/>'
            f'<ellipse cx="27" cy="66" rx="7" ry="9" transform="rotate(-18 27 66)"'
            f' fill="#ffffff" opacity="0.85"/>'
            f'<path d="M44,20 L54,16" stroke="#ffffff" stroke-width="2.6"'
            f' opacity="0.8" stroke-linecap="round"/>'
            f'<path d="M41,27 L51,23" stroke="#ffffff" stroke-width="2.6"'
            f' opacity="0.8" stroke-linecap="round"/>')


def icon_trypanosome(c):
    centre = lambda t: (8 + 84 * t, 56 + 16 * math.sin(2 * math.pi * 0.9 * t + 0.6))
    body = ribbon(centre, taper(peak=8.5, head=0.30, tail=0.05))
    # free flagellum running off the anterior end
    flag = ("M76,44 C86,36 92,30 96,22")
    # undulating membrane along the body
    memb = ("M16,52 C30,34 44,60 58,40 C68,26 76,32 80,40")
    return (f'<path d="{body}" fill="{c}"/>'
            f'<path d="{memb}" stroke="{c}" stroke-width="3" fill="none"'
            f' opacity="0.55" stroke-linecap="round"/>'
            f'<path d="{flag}" stroke="{c}" stroke-width="3.4" fill="none"'
            f' stroke-linecap="round"/>'
            f'<ellipse cx="44" cy="58" rx="8" ry="6.5" fill="#ffffff" opacity="0.85"/>'
            f'<circle cx="66" cy="50" r="3.2" fill="#ffffff" opacity="0.85"/>')


def icon_amoeba(c):
    body = ("M28,20 C44,8 62,12 70,24 C78,36 94,34 92,48 C90,62 76,62 74,74 "
            "C72,88 54,96 40,88 C26,80 30,68 20,62 C8,54 12,30 28,20 Z")
    return (f'<path d="{body}" fill="{c}"/>'
            f'<circle cx="46" cy="50" r="11" fill="#ffffff" opacity="0.85"/>'
            f'<circle cx="46" cy="50" r="4" fill="{c}" opacity="0.55"/>'
            f'<circle cx="68" cy="66" r="5" fill="#ffffff" opacity="0.5"/>')


def icon_spore(c):
    coil = []
    for i in range(5):                       # polar tube coiled inside the spore
        y = 40 + i * 9
        coil.append(f'<path d="M32,{y} C44,{y - 7} 58,{y - 7} 68,{y}"'
                    f' stroke="#ffffff" stroke-width="3" fill="none" opacity="0.8"'
                    f' stroke-linecap="round"/>')
    return (f'<ellipse cx="50" cy="50" rx="27" ry="42" fill="{c}"/>'
            f'<ellipse cx="50" cy="28" rx="9" ry="7" fill="#ffffff" opacity="0.55"/>'
            + "".join(coil))


def icon_tick(c):
    legs = []
    for y, knee, foot in [(40, 12, 4), (50, 15, 16), (60, 15, 28), (70, 12, 38)]:
        legs.append(f'<path d="M32,{y} L{32 - knee},{y - 6} L{32 - knee - 10},{foot}"'
                    f' stroke="{c}" stroke-width="3.4" fill="none"'
                    f' stroke-linecap="round" stroke-linejoin="round"/>')
        legs.append(f'<path d="M68,{y} L{68 + knee},{y - 6} L{68 + knee + 10},{foot}"'
                    f' stroke="{c}" stroke-width="3.4" fill="none"'
                    f' stroke-linecap="round" stroke-linejoin="round"/>')
    return ("".join(legs)
            + f'<ellipse cx="50" cy="60" rx="26" ry="31" fill="{c}"/>'
            + f'<ellipse cx="50" cy="60" rx="15" ry="20" fill="#ffffff" opacity="0.22"/>'
            + f'<ellipse cx="50" cy="30" rx="12" ry="10" fill="{c}"/>'
            + f'<path d="M44,21 L41,10 M56,21 L59,10 M50,20 L50,7" stroke="{c}"'
              f' stroke-width="3.2" stroke-linecap="round"/>')


ICONS = {"nematode": icon_nematode, "fluke": icon_fluke, "tapeworm": icon_tapeworm,
         "apicomplexan": icon_apicomplexan, "trypanosome": icon_trypanosome,
         "amoeba": icon_amoeba, "spore": icon_spore, "tick": icon_tick}


# --- host silhouettes (off by default; --host-icons brings them back) --------
# Real silhouettes from PhyloPic (https://www.phylopic.org), all released under
# CC0, cached in images/silhouettes/ so the script needs no network.
#   human  f36c3c31-8b66-4413-a39f-a2b2b422e7ff  T. Michael Keesey
#   rat    6b8ecf3f-a5c2-4ca4-adbc-1c7280c380d4  Arcadia Science
#   mouse  36dc0476-ae7d-49ed-85c4-220139930bfc  Cagri Cevrim
#   pig    216ee85d-1696-49ae-a62c-d1da6fef06fc  anonymous
SILHOUETTE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "images", "silhouettes")
SILHOUETTE_CREDIT = ("Host silhouettes: PhyloPic (CC0) — T. M. Keesey, "
                     "Arcadia Science, C. Cevrim")
# file, mirror so every animal faces right, and size relative to the others
HOST_SILHOUETTE = {
    9606:  dict(file="human.svg", flip=True,  rel=1.00),
    10116: dict(file="rat.svg",   flip=False, rel=1.00),
    10090: dict(file="mouse.svg", flip=False, rel=0.72),
    9823:  dict(file="pig.svg",   flip=True,  rel=1.00),
}


@functools.lru_cache(maxsize=None)
def load_silhouette(filename):
    """Return (viewBox, drawing) for a cached PhyloPic SVG, with fills stripped."""
    with open(os.path.join(SILHOUETTE_DIR, filename)) as fh:
        svg = fh.read()
    vb = [float(v) for v in re.search(r'viewBox="([^"]+)"', svg).group(1).replace(",", " ").split()]
    body = svg[svg.index(">", svg.index("<svg")) + 1: svg.rindex("</svg>")]
    body = re.sub(r"<metadata>.*?</metadata>", "", body, flags=re.S)
    body = re.sub(r'\s(fill|stroke)="[^"]*"', "", body)
    return tuple(vb), body.strip()


def place_silhouette(taxid, color, x, y, box_w, box_h):
    """Fit a host silhouette into (x, y, box_w, box_h), bottom-aligned."""
    spec = HOST_SILHOUETTE[taxid]
    (vx, vy, vw, vh), body = load_silhouette(spec["file"])
    k = min(box_w / vw, box_h / vh) * spec["rel"]
    dw, dh = vw * k, vh * k
    ox = x + (box_w - dw) / 2.0
    oy = y + (box_h - dh)                       # sit on a shared baseline
    if spec["flip"]:
        t = (f"translate({ox + dw:.2f},{oy:.2f}) scale({-k:.5f},{k:.5f}) "
             f"translate({-vx:.2f},{-vy:.2f})")
    else:
        t = (f"translate({ox:.2f},{oy:.2f}) scale({k:.5f}) "
             f"translate({-vx:.2f},{-vy:.2f})")
    return f'<g transform="{t}" fill="{color}">{body}</g>'


# --- text --------------------------------------------------------------------
def text(x, y, s, size=16, color=INK, weight="normal", style="normal",
         anchor="start", spacing="0"):
    s = (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}"'
            f' fill="{color}" font-weight="{weight}" font-style="{style}"'
            f' text-anchor="{anchor}" letter-spacing="{spacing}">{s}</text>')


def place_icon(draw, color, x, y, size):
    """Place a 100x100 icon with its top-left corner at (x, y)."""
    k = size / 100.0
    return (f'<g transform="translate({x:.1f},{y:.1f}) scale({k:.4f})">'
            f'{draw(color)}</g>')


# --- the figure --------------------------------------------------------------
def build(cfg, title, host_icons=False):
    hosts = cfg["hosts"]
    parasites = cfg["parasites"]
    gcolors = cfg["parasite_groups"]

    by_group = {g: [] for g in GROUP_ORDER}
    for taxid, p in parasites.items():
        by_group.setdefault(p["group"], []).append((p["label"], set(p["hosts"])))
    for g in by_group:
        by_group[g].sort()

    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}"'
           f' viewBox="0 0 {W} {H}">',
           f'<rect width="{W}" height="{H}" fill="#ffffff"/>']

    # --- title
    y = MARGIN + 8
    out.append(text(MARGIN, y, title, size=30, weight="600"))
    out.append(text(MARGIN, y + 28,
                    f"{len(parasites)} parasite species in {len(hosts)} host species",
                    size=17, color=MUTED))

    # --- host band, across the full width
    band_y = y + 74
    gap = 20
    card_w = (W - 2 * MARGIN - 3 * gap) / 4.0
    card_h = 138 if host_icons else 104
    out.append(text(MARGIN, band_y - 12, "HOSTS", size=13, color=MUTED,
                    weight="600", spacing="1.6"))
    for i, taxid in enumerate(HOST_ORDER):
        h = hosts[taxid]
        c = h["color"]
        x = MARGIN + i * (card_w + gap)
        n = sum(1 for p in parasites.values() if taxid in p["hosts"])
        out.append(f'<rect x="{x:.1f}" y="{band_y}" width="{card_w:.1f}"'
                   f' height="{card_h}" rx="14" fill="{c}" opacity="0.08"/>')
        out.append(f'<rect x="{x:.1f}" y="{band_y + 14}" width="5"'
                   f' height="{card_h - 28}" rx="2.5" fill="{c}"/>')
        if host_icons:
            out.append(place_silhouette(taxid, c, x + 16, band_y + 12, 130, 114))
            tx, ty = x + 156, band_y + 52
        else:
            tx, ty = x + 28, band_y + 42
        out.append(text(tx, ty, COMMON_NAME[taxid], size=22, weight="600"))
        out.append(text(tx, ty + 24, h["label"], size=14.5, color=MUTED,
                        style="italic"))
        out.append(text(tx, ty + 56, f"{n} parasites", size=15, color=c,
                        weight="600"))

    # --- parasite panels
    top = band_y + card_h + 56
    out.append(text(MARGIN, top - 16, "PARASITES", size=13,
                    color=MUTED, weight="600", spacing="1.6"))
    col_gap = 22
    unit = (W - 2 * MARGIN - 4 * col_gap) / 5.0
    columns = [
        (2, ["Nematoda"]),
        (1, ["Trematoda", "Cestoda"]),
        (1, ["Apicomplexa", "Kinetoplastida"]),
        (1, ["Other protozoa", "Microsporidia"]),
    ]
    head_h, line_h, pad_b = 84, 30, 16

    x = MARGIN
    for span, groups in columns:
        w = span * unit + (span - 1) * col_gap
        yy = top
        for g in groups:
            members = by_group[g]
            ncol = 2 if span == 2 else 1
            rows = math.ceil(len(members) / ncol)
            ph = head_h + rows * line_h + pad_b
            c = gcolors[g]
            out.append(f'<rect x="{x}" y="{yy}" width="{w:.1f}" height="{ph}" rx="14"'
                       f' fill="{PANEL_BG}" stroke="{PANEL_LINE}"/>')
            out.append(f'<rect x="{x}" y="{yy}" width="{w:.1f}" height="5" rx="2.5"'
                       f' fill="{c}"/>')
            out.append(place_icon(ICONS[GROUP_ICON[g]], c, x + 14, yy + 16, 46))
            out.append(text(x + 72, yy + 40, g, size=19, weight="600"))
            out.append(text(x + 72, yy + 58, GROUP_SUB[g], size=12.5, color=MUTED))
            out.append(text(x + w - 16, yy + 44, str(len(members)), size=23,
                            color=c, weight="700", anchor="end"))

            colw = (w - 24) / ncol
            for i, (name, phosts) in enumerate(members):
                ci, ri = (i // rows, i % rows) if ncol == 2 else (0, i)
                nx = x + 14 + ci * colw
                ny = yy + head_h + ri * line_h
                for j, (htax, _) in enumerate(DOT_HOSTS):
                    dx = nx + 5 + j * 11
                    if htax in phosts:
                        out.append(f'<circle cx="{dx}" cy="{ny - 4.5}" r="4.6"'
                                   f' fill="{hosts[htax]["color"]}"/>')
                out.append(text(nx + 54, ny, name, size=14.5, style="italic"))
            yy += ph + 20
        x += w + col_gap

    # --- legend
    ly = H - 30
    out.append(text(MARGIN, ly, "Host:", size=13, color=MUTED, weight="600"))
    lx = MARGIN + 44
    for htax, label in DOT_HOSTS:
        out.append(f'<circle cx="{lx}" cy="{ly - 4.5}" r="4.6"'
                   f' fill="{hosts[htax]["color"]}"/>')
        out.append(text(lx + 12, ly, label, size=13, color=MUTED))
        lx += 30 + 7.2 * len(label)
    if host_icons:
        out.append(text(W - MARGIN, ly, SILHOUETTE_CREDIT, size=11.5, color=MUTED,
                        anchor="end"))
    out.append("</svg>")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("-c", "--config", default=os.path.join(root, "config.yml"))
    ap.add_argument("-o", "--output", default=os.path.join(root, "docs", "hosts_parasites.svg"))
    ap.add_argument("-t", "--title", default="Hosts and parasites in OrthoHPI 2.0")
    ap.add_argument("--host-icons", action="store_true",
                    help="draw the PhyloPic host silhouettes in the host cards")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    svg = build(cfg, args.title, host_icons=args.host_icons)
    with open(args.output, "w") as fh:
        fh.write(svg)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
