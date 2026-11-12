#!/usr/bin/env python3
"""Figure 1, the pbzarr data model, drawn flat in the manner of AnnData's Figure 1.

A store holds one reference and any number of tracks. Every track runs along
the same position axis, every base of every contig end to end, so the tracks
line up row for row. contig and contig_length hold the reference: one entry
per contig, drawn here across the rows of its contig. A track has the position
axis and at most one labeled axis (sample, or any other column), whose labels
sit in a small array above it.

Each track is cut into chunks along position (solid lines), so contig edges
(dashed) can fall inside a chunk. One chunk (red) is
encoded by the array's codec pipeline and stored under its key in the store on
the right.

Writes figures/figure1.png and .svg, drawn at print size to fit a 160 mm text
width (6-7 pt text), so the SVG is placed in Illustrator at 100%. Earlier drafts are in
figures/archive/ (figure1-vcfzarr-style, figure1-sheets). Run from the repo root:
    pixi run python figures/figure1.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.transforms import Bbox

ROOT = Path(__file__).resolve().parent.parent
INK = "#444"
MONO = "Courier New"
PBZ = "#C0304F"  # PBZ's color in figures 2 and 3: here, the chunk we follow
# Overview colors: the reference in blue, labels in orange, one hue per track.
REFERENCE, LABELS = "#B5D4F4", "#F5C4B3"
CELL = 0.15  # one position row, inches
COL = 0.24  # one track column, inches: wide enough for a value
ROWS = 9  # position rows, the last two (the contigs after chr2) drawn as an ellipsis
SHOWN = 7  # rows drawn as cells
CHUNK = 3  # rows per chunk
CONTIGS = [("chr1", "248,956,422", 0, 4), ("chr2", "242,193,529", 4, 7), ("…", "…", 7, 9)]
WIDTH, HEIGHT = 180 / 25.4, 5.0

plt.rcParams.update({"font.family": "Arial", "font.size": 7, "svg.fonttype": "none"})


def block(ax, x, top, cols, rows, color, cell_w=CELL):
    """Flat cells with a darker outline; top is the upper edge."""
    w, h = cols * cell_w, rows * CELL
    ax.add_patch(Rectangle((x, top - h), w, h, fc=color, ec="none", zorder=2))
    for c in range(1, cols):
        ax.plot([x + c * cell_w] * 2, [top - h, top], color="white", lw=0.5, zorder=3)
    for r in range(1, rows):
        ax.plot([x, x + w], [top - r * CELL] * 2, color="white", lw=0.5, zorder=3)
    ax.add_patch(Rectangle((x, top - h), w, h, fc="none", ec=INK, lw=0.6, zorder=4))


def name(ax, x, y, s, va="baseline"):
    ax.text(x, y, s, ha="center", va=va, fontsize=6.5, family=MONO, weight="bold")


def arrow(ax, start, end, color=INK, lw=0.6):
    ax.annotate("", end, start, zorder=4, arrowprops=dict(
        arrowstyle="-|>", mutation_scale=5 if lw < 1 else 7, shrinkA=0, shrinkB=0, lw=lw, color=color))


def dots(ax, x, y):
    """A vertical ellipsis as three dots (Arial has none)."""
    ax.plot([x] * 3, [y - 0.035, y, y + 0.035], ls="none", marker="o", ms=0.6, color=INK, zorder=5)


def frame(ax, x, y, w, h):
    """A dotted frame with its lower left corner at x, y."""
    ax.add_patch(Rectangle((x, y), w, h, fc="none", ec=INK, lw=0.6, ls=(0, (1, 1.5)), zorder=1))


fig = plt.figure(figsize=(WIDTH, HEIGHT))
ax = fig.add_axes([0, 0, 1, 1])
ax.set(xlim=(0, WIDTH), ylim=(0, HEIGHT), aspect="equal")
ax.axis("off")

TOP = HEIGHT - 1.0
BOTTOM = TOP - ROWS * CELL
NAMES = TOP + 0.25 + CELL  # the row of array and track names

# The reference: one entry per contig, spanning that contig's rows.
x = 0.1
for column, RW in [("contig", 0.36), ("contig_length", 0.64)]:
    name(ax, x + RW / 2, NAMES, column)
    for contig, length, start, end in CONTIGS:
        y0, y1 = TOP - end * CELL, TOP - start * CELL
        if contig == "…":
            dots(ax, x + RW / 2, (y0 + y1) / 2)
            continue
        ax.add_patch(Rectangle((x, y0), RW, y1 - y0, fc=REFERENCE, ec=INK, lw=0.6, zorder=2))
        ax.text(x + RW / 2, (y0 + y1) / 2, contig if column == "contig" else length, ha="center", va="center",
                fontsize=6, family=MONO, zorder=5)
    x += RW + 0.08
ref_right = x - 0.08

# The position axis.
axis_x = ref_right + 0.1
arrow(ax, (axis_x, TOP), (axis_x, BOTTOM))
ax.text(axis_x, BOTTOM - 0.05, "position", ha="center", va="top", fontsize=7)

# The tracks: (name, color, labeled axis, its labels, and the values of the
# rows drawn). Coverage is an integer, DNase a float, repeats a boolean. The
# last one holds the chunk we follow, so its arrow crosses no other track.
tracks = [
    ("track_cov", "#CECBF6", "sample", ["s1", "s2", "s3"],
     [[31, 28, 35], [33, 30, 34], [30, 27, 36], [12, 9, 14], [29, 31, 33], [34, 29, 30], [32, 30, 31]]),
    ("track_dnase", "#9FE1CB", "sample", ["s1", "s2", "s3"],
     [[0.2, 0.1, 0.3], [1.8, 2.1, 1.5], [2.4, 2.6, 2.2], [0.1, 0.0, 0.2], [0.5, 0.4, 0.7], [0.3, 0.2, 0.1],
      [1.2, 0.9, 1.4]]),
    ("track_repeats", "#F4C0D1", "repeat", ["LINE", "Alu"],
     [[0, 0], [0, 1], [0, 1], [1, 0], [1, 0], [0, 0], [0, 1]]),
]


def draw_track(x, top, names_y, track, color, dim, entries, values):
    """A track at x with its cells' top at top: name, labeled axis, values and
    the ellipses for the rows not drawn. Returns its width."""
    cols = len(entries)
    w = cols * COL
    block(ax, x, top, cols, SHOWN, color, cell_w=COL)
    for r, row in enumerate(values):
        for c, v in enumerate(row):
            text = "TF"[1 - v] if track == "track_repeats" else f"{v:.1f}" if track == "track_dnase" else f"{v}"
            ax.text(x + (c + 0.5) * COL, top - (r + 0.5) * CELL, text, ha="center", va="center", fontsize=5.5,
                    family=MONO, zorder=5)
    for c in range(cols):  # the rows after chr2
        dots(ax, x + (c + 0.5) * COL, top - (SHOWN + 1) * CELL)
    name(ax, x + w / 2, names_y, track)
    block(ax, x, top + 0.08 + CELL, cols, 1, LABELS, cell_w=COL)
    for c, entry in enumerate(entries):
        ax.text(x + (c + 0.5) * COL, top + 0.08 + CELL / 2, entry, ha="center", va="center", fontsize=6,
                family=MONO, zorder=5)
    ax.text(x + w / 2, top + 0.11 + CELL, dim, ha="center", va="bottom", fontsize=6.5, family=MONO)
    return w


x = axis_x + 0.16
for track, *rest in tracks:
    w = draw_track(x, TOP, NAMES, track, *rest)
    for r in range(CHUNK, SHOWN, CHUNK):  # chunk edges
        ax.plot([x, x + w], [TOP - r * CELL] * 2, color=INK, lw=0.6, zorder=4)
    x += w + 0.26
right = x - 0.26
chunk_w = len(tracks[-1][3]) * COL
chunk_x = right - chunk_w

# Contig edges run straight across the reference and every track.
for _, _, start, _ in CONTIGS[1:]:
    y = TOP - start * CELL
    ax.plot([0.04, right + 0.08], [y, y], color=INK, lw=0.5, ls=(0, (3, 2)), zorder=6)

# The store bracket, named like the store's root in the tree.
by = NAMES + 0.17  # above the array names
ax.plot([0.04, 0.04, right + 0.08, right + 0.08], [by - 0.05, by, by, by - 0.05], color=INK, lw=0.6)
ax.text((0.04 + right + 0.08) / 2, by + 0.05, "store.pbz", ha="center", va="bottom", fontsize=7, family=MONO)

# B, below A and on its own: tracks that share position combine in one array
# operation. track_dnase times track_repeats keeps each sample's DNase only
# where a repeat is (T), once per repeat, so the result runs over position,
# repeat and sample.
values = {t: v for t, _, _, _, v in tracks}
repeats, dnase = values["track_repeats"], values["track_dnase"]
B_NAMES = BOTTOM - 0.45
B_TOP = B_NAMES - 0.45  # top of B's cells
SIGN, GAP_B = 0.3, 0.12  # room for the x and = signs; gap between the repeats' blocks
b_width = 3 * COL + SIGN + 2 * COL + SIGN + 6 * COL + GAP_B
x = 0.04 + (right + 0.04 - b_width) / 2  # centered under A
for track, sign in [("track_dnase", "×"), ("track_repeats", "=")]:
    x += draw_track(x, B_TOP, B_NAMES, *next(t for t in tracks if t[0] == track))
    ax.text(x + SIGN / 2, B_TOP - SHOWN * CELL / 2, sign, ha="center", va="center", fontsize=10)
    x += SIGN
# The result as one 3D block: position down, sample across, repeat into the
# page. The Alu slice sits behind the LINE slice, offset up and to the right.
DX, DY = 0.32, 0.26  # offset of each slice behind the one in front
name(ax, x + 1.5 * COL + DX / 2, B_NAMES + DY / 2, "dnase_by_repeat")
for k in (1, 0):  # back slice first; each slice's layers sit above the one behind
    repeat = ["LINE", "Alu"][k]
    bx, top = x + k * DX, B_TOP + k * DY
    z = 3 - 1.5 * k
    ax.add_patch(Rectangle((bx, top - SHOWN * CELL), 3 * COL, SHOWN * CELL, fc="#eeeeee", ec="none", zorder=z))
    for r in range(SHOWN):
        y = top - (r + 0.5) * CELL
        if repeats[r][k]:
            ax.add_patch(Rectangle((bx, y - CELL / 2), 3 * COL, CELL, fc="#9FE1CB", ec="none", zorder=z + 0.1))
            for j, v in enumerate(dnase[r]):
                ax.text(bx + (j + 0.5) * COL, y, f"{v:.1f}", ha="center", va="center", fontsize=5.5, family=MONO,
                        zorder=z + 0.3)
    for c in range(1, 3):
        ax.plot([bx + c * COL] * 2, [top - SHOWN * CELL, top], color="white", lw=0.5, zorder=z + 0.2)
    for r in range(1, SHOWN):
        ax.plot([bx, bx + 3 * COL], [top - r * CELL] * 2, color="white", lw=0.5, zorder=z + 0.2)
    ax.add_patch(Rectangle((bx, top - SHOWN * CELL), 3 * COL, SHOWN * CELL, fc="none", ec=INK, lw=0.6,
                           zorder=z + 0.4))
    # The repeat's label, past the slice's bottom right corner.
    ax.text(bx + 3 * COL + 0.05, top - SHOWN * CELL - 0.02, repeat, ha="left", va="top", fontsize=6, family=MONO,
            bbox=dict(boxstyle="square,pad=0.15", fc=LABELS, ec=INK, lw=0.6), zorder=6)
# The block's edges, joining the slices' corners.
for cx, cy in [(x, B_TOP), (x + 3 * COL, B_TOP), (x + 3 * COL, B_TOP - SHOWN * CELL)]:
    ax.plot([cx, cx + DX], [cy, cy + DY], color=INK, lw=0.6, zorder=1)
# The sample labels, above the front slice, and the ellipses below it.
block(ax, x, B_TOP + 0.08 + CELL, 3, 1, LABELS, cell_w=COL)
for j, sample in enumerate(["s1", "s2", "s3"]):
    ax.text(x + (j + 0.5) * COL, B_TOP + 0.08 + CELL / 2, sample, ha="center", va="center", fontsize=6,
            family=MONO, zorder=5)
    dots(ax, x + (j + 0.5) * COL, B_TOP - (SHOWN + 1) * CELL)
ax.text(x + 3 * COL + DX / 2 + 0.25, B_TOP - SHOWN * CELL + DY / 2 - 0.08, "repeat", ha="left", va="center", fontsize=6.5, family=MONO,
        rotation=39, color=INK)

# C, to the right of A: the chunk we follow, track_repeats' second chunk.
c_top, c_bottom = TOP - CHUNK * CELL, TOP - 2 * CHUNK * CELL
ax.add_patch(Rectangle((chunk_x, c_bottom), chunk_w, c_top - c_bottom, fc=PBZ, alpha=0.25, ec="none", zorder=3))
ax.add_patch(Rectangle((chunk_x, c_bottom), chunk_w, c_top - c_bottom, fc="none", ec=PBZ, lw=1.4, zorder=7))
mid = (c_top + c_bottom) / 2

# It goes through the array's codec pipeline, as in Zarr v3: optional
# array -> array codecs (filters), one array -> bytes codec, then optional
# bytes -> bytes codecs (compressors). The codecs are chosen per array.
CODECS = [("Filter (array → array)", "e.g. transpose, delta"),
          ("Serialize (array → bytes)", "e.g. raw bytes, sharding"),
          ("Compress (bytes → bytes)", "e.g. Zstandard, Blosc")]
COMP_X, COMP_W, BOX_W, BOX_H, GAP = right + 0.3, 1.5, 1.2, 0.27, 0.4
box_x = COMP_X + (COMP_W - BOX_W) / 2
box_ys = [TOP + 0.15 - k * GAP for k in range(len(CODECS))]  # box centers, the first level with the label arrays
codec_bottom = box_ys[-1] - BOX_H / 2 - 0.1
frame(ax, COMP_X, codec_bottom, COMP_W, box_ys[0] - box_ys[-1] + BOX_H + 0.32)
ax.text(COMP_X + COMP_W / 2, box_ys[0] + BOX_H / 2 + 0.1, "Codec pipeline", ha="center", va="center", fontsize=7)
for (kind, example), by_ in zip(CODECS, box_ys):
    ax.add_patch(FancyBboxPatch((box_x, by_ - BOX_H / 2), BOX_W, BOX_H, boxstyle="round,pad=0,rounding_size=0.04",
                                fc="white", ec=INK, lw=0.6, zorder=3))
    ax.text(box_x + BOX_W / 2, by_ + 0.045, kind, ha="center", va="center", fontsize=6.5, zorder=5)
    ax.text(box_x + BOX_W / 2, by_ - 0.06, example, ha="center", va="center", fontsize=6, color="#222", zorder=5)
for upper, lower in zip(box_ys, box_ys[1:]):
    arrow(ax, (box_x + BOX_W / 2, upper - BOX_H / 2), (box_x + BOX_W / 2, lower + BOX_H / 2), color=PBZ, lw=1.2)
# From the chunk out, up and into the first codec, with right-angle bends.
bend_x = (chunk_x + chunk_w + COMP_X) / 2
ax.plot([chunk_x + chunk_w, bend_x, bend_x], [mid, mid, box_ys[0]], color=PBZ, lw=1.2, zorder=4,
        solid_joinstyle="miter")
arrow(ax, (bend_x, box_ys[0]), (box_x, box_ys[0]), color=PBZ, lw=1.2)
out_y = box_ys[-1]

# ...and into storage, stacked below, under its key: the chunk's index along
# position, then along the track's labeled axis.
# Format 0.5 (pbzarr-rs docs/STORE_SPEC.md): every array sits at the root, with
# no track groups; labels are shared by dimension name.
tree = [(0, "store.pbz/"), (1, "contig/"), (1, "contig_length/"), (1, "sample/"), (1, "repeat/"),
        (1, "track_cov/"), (1, "track_dnase/"), (1, "track_repeats/"), (2, "c/0/0"), (2, "c/1/0"), (2, "c/2/0")]
KEY = "c/1/0"
LINE = 0.13
STEP = 0.18  # indent per level
top = codec_bottom - 0.15 - 0.22  # 0.22: the frame's room above the first entry
frame(ax, COMP_X, top - len(tree) * LINE - 0.08, COMP_W, len(tree) * LINE + 0.3)
ax.text(COMP_X + COMP_W / 2, top + 0.11, "File system or cloud storage", ha="center", va="center", fontsize=7)
ys = {}
for k, (level, entry) in enumerate(tree):
    tx, ty = COMP_X + 0.15 + level * STEP, top - (k + 0.6) * LINE
    ys[k] = ty
    if level:  # a branch from the parent's line
        ax.plot([tx - 0.11, tx - 0.03], [ty, ty], color="black", lw=0.6, zorder=2)
    red = entry == KEY
    ax.text(tx, ty, entry, ha="left", va="center", fontsize=6.5, family=MONO, zorder=5,
            color=PBZ if red else "black", weight="bold")
for k, (level, _) in enumerate(tree):  # vertical lines, from a parent to its last child
    children = [j for j in range(k + 1, len(tree)) if tree[j][0] == level + 1]
    children = [j for j in children if all(tree[i][0] > level for i in range(k + 1, j + 1))]
    if children:
        lx = COMP_X + 0.15 + (level + 1) * STEP - 0.11
        ax.plot([lx, lx], [ys[k] - 0.05, ys[children[-1]]], color="black", lw=0.6, zorder=2)
key_y = ys[[e for _, e in tree].index(KEY)]
# From the last codec out to the right, down, and back left to the key.
turn_x = COMP_X + COMP_W - 0.07
ax.plot([box_x + BOX_W, turn_x, turn_x], [out_y, out_y, key_y], color=PBZ, lw=1.2, zorder=4, solid_joinstyle="miter")
arrow(ax, (turn_x, key_y), (COMP_X + 0.15 + 2 * STEP + 0.33, key_y), color=PBZ, lw=1.2)

# Panel letters.
for letter, lx, ly in [("A", 0.02, NAMES + 0.42), ("B", 0.02, B_NAMES + 0.2), ("C", COMP_X - 0.05, NAMES + 0.42)]:
    ax.text(lx, ly, letter, ha="left", va="bottom", fontsize=10, weight="bold")

# Crop to what is drawn: the axes span the whole canvas, so "tight" alone keeps its empty edges.
renderer = fig.canvas.get_renderer()
drawn = Bbox.union([a.get_window_extent(renderer) for a in ax.patches + ax.lines + ax.texts])
crop = drawn.transformed(fig.dpi_scale_trans.inverted()).padded(0.04)
print(f"Figure 1 is {crop.width * 25.4:.0f} x {crop.height * 25.4:.0f} mm")
fig.savefig(ROOT / "figures" / "figure1.png", dpi=600, bbox_inches=crop)
fig.savefig(ROOT / "figures" / "figure1.svg", bbox_inches=crop)
