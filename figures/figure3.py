#!/usr/bin/env python3
"""Figure 3, scans: mean depth from the cohort files of figure 2.

A: how each format computes mean depth, over the whole genome or (with a BED
   file) over regions.
B: short-read sample scaling of the whole-genome and region scans.
C, D: whole-genome and region scan wall time of the largest cohort of each
   data type, one panel per data type.

The over10 scans and the scan memory are left out. Reads the scans/*.tsv
files copied from Redwood into results/ (3 cold-cache trials per scan; the
figure uses their median) and writes figures/figure3.png and .svg (the whole
figure), plus each panel as an SVG for Illustrator in figures/panels/
(figure3-A to figure3-D, legend). The figure is drawn at print size for an A4
page (180 mm wide, 6-7 pt text), so it and the panels are placed at 100%. Run
from the repo root:
    pixi run python figures/figure3.py
"""

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Patch
from matplotlib.transforms import Bbox, blended_transform_factory

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
PANELS = ROOT / "figures" / "panels"
# As in figure 2: BigWig, D4, PBZ, with PBZ the vibrant color.
TOOLS = {"bigwig": ("BigWig", "#9CC9C0"), "d4": ("D4", "#A7BEDB"), "pbz": ("PBZ", "#C0304F")}
DARK = {"bigwig": "#3F8073", "d4": "#4A6E9C", "pbz": "#C0304F"}
INK = "#444"
MONO = "Courier New"
SCANS = {"full_scan": "Calculating mean coverage of the whole genome",
         "region_scan": "Calculating mean coverage of 10,000 random regions"}


def read(pattern):
    """One row per scan: the median of its trials, with N, tool and scan from the path."""
    rows = []
    for f in sorted(RESULTS.glob(pattern)):
        tool, scan = f.stem.split(".")
        trials = pd.read_csv(f, sep="\t")
        rows.append({"n": int(re.search(r"/n(\d+)/", str(f)).group(1)), "tool": tool, "scan": scan,
                     "s": trials.s.median()})
    return pd.DataFrame(rows)


# Short read: PBZ and D4/BigWig ran as separate jobs; only the first build was scanned.
sr = pd.concat([read("cohort-pbz/n*/b1/scans/*_scan.tsv"), read("cohort-d4bw/n*/b1/scans/*_scan.tsv")])

largest = pd.concat([
    sr[sr.n == 1024].assign(dataset="Short read"),
    read("cohort-ont/n1019/scans/*_scan.tsv").assign(dataset="Long read"),
    read("cohort-rna/n60/scans/*_scan.tsv").assign(dataset="Long-read RNA"),
])
DATASETS = list(largest.dataset.unique())


def seconds_text(s):
    return f"{s:,.0f} s" if s >= 10 else f"{s:.1f} s"


def letter(ax, text, dx=-22, dy=6):
    return ax.annotate(text, (0, 1), xycoords="axes fraction", xytext=(dx, dy), textcoords="offset points",
                       ha="right", va="bottom", fontsize=10, weight="bold")


def save_panels(fig, panels, letters):
    """Save each panel (its axes, or the legend) as its own SVG, cut out of fig
    with everything else and the panel letters hidden. Text stays text."""
    renderer = fig.canvas.get_renderer()
    for artist in letters:
        artist.set_visible(False)
    for name, parts in panels.items():
        for artist in fig.axes + fig.legends:
            artist.set_visible(artist in parts)
        box = Bbox.union([p.get_tightbbox(renderer) for p in parts])
        fig.savefig(PANELS / f"{name}.svg", bbox_inches=box.transformed(fig.dpi_scale_trans.inverted()).padded(0.05))


# The schematic is drawn in units of UNIT inches, so its boxes and text keep
# their size whatever the size of the panel.
UNIT = 0.42


def box(ax, x, y, text, tool=False, w=1.1, h=0.38, dashed=False):
    """A file (white) or a tool (grey, monospace) centered on x, y."""
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0,rounding_size=0.08",
                                fc="#eeeeee" if tool else "white", ec=INK, lw=0.6, ls="--" if dashed else "-",
                                zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=5.5 if tool else 6, family=MONO if tool else None,
            zorder=4, linespacing=1)


def arrow(ax, start, end):
    ax.annotate("", end, start, zorder=2, arrowprops=dict(
        arrowstyle="-|>", mutation_scale=5, shrinkA=0, shrinkB=0, lw=0.6, color=INK))


def schematic(ax, width, height):
    """One column per tool: the cohort file(s), the tool that reads them, and
    the mean depth table(s) it writes. A BED file of regions, used only by the
    region scan, feeds every tool along the dashed line."""
    w, h = width / UNIT, height / UNIT
    ax.set(xlim=(0, w), ylim=(0, h), aspect="equal")
    ax.axis("off")
    head, data, scan, out = h - 0.3, h - 0.95, h - 1.8, h - 2.65
    for y, label in [(data, "Input"), (out, "Output")]:
        ax.text(0.05, y, label, va="center", fontsize=7, color=INK)
    box(ax, 0.85, scan, "regions.bed\n(region scan)", w=1.45, h=0.55, dashed=True)
    col = (w - 1.7) / 3
    centers = [1.7 + col * (i + 0.5) for i in range(3)]
    ax.plot([1.58, centers[-1]], [scan, scan], color="#999", lw=0.6, ls="--", zorder=1)
    for cx, (tool, (name, color)) in zip(centers, TOOLS.items()):
        ax.add_patch(FancyBboxPatch((cx - col / 2 + 0.06, 0), col - 0.12, h, boxstyle="round,pad=0,rounding_size=0.12",
                                    fc=color, alpha=0.25, lw=0, zorder=0))
        ax.text(cx, head, name, ha="center", va="center", fontsize=7, weight="bold", color=DARK[tool])
        if tool == "bigwig":
            for x, s in [(cx - 1.0, "s1"), (cx + 1.0, "sN")]:
                box(ax, x, data, f"{s}.bw")
                arrow(ax, (x, data - 0.19), (x, scan + 0.19))
                box(ax, x, scan, "pyBigWig", tool=True, w=1.28)
                arrow(ax, (x, scan - 0.19), (x, out + 0.19))
                box(ax, x, out, f"{s}.tsv")
            for y in (data, out):
                ax.text(cx, y, "…", ha="center", va="center", fontsize=6)
            continue
        command = {"d4": "d4tools stat", "pbz": "pbz stat"}[tool]
        box(ax, cx, data, f"cohort.{tool}", w=1.3)
        arrow(ax, (cx, data - 0.19), (cx, scan + 0.19))
        box(ax, cx, scan, command, tool=True, w=1.6)
        arrow(ax, (cx, scan - 0.19), (cx, out + 0.19))
        box(ax, cx, out, "means.tsv", w=1.3)


def scaling(ax, scan):
    """One line per tool over cohort size, with each tool's time at N = 1,024
    written past the line's end, the fastest in bold."""
    rows = sr[sr.scan == scan]
    for tool, (_, color) in TOOLS.items():
        line = rows[rows.tool == tool].sort_values("n")
        ax.plot(line.n, line.s, color=color, lw=1.6 if tool == "pbz" else 1, marker="o", ms=2.5, zorder=3)
    ax.set_title(SCANS[scan], loc="left", weight="bold")
    ax.set_ylabel("Wall time (s)")
    ax.set_xlabel("Samples")
    ax.set_xticks([0, 256, 512, 768, 1024])
    ax.set(xlim=(0, 1400), ylim=(0, None))

    # Spread the end labels so they are at least a line apart.
    ends = rows[rows.n == 1024].set_index("tool").s.sort_values()
    to_axes = ax.transData + ax.transAxes.inverted()
    at = blended_transform_factory(ax.transData, ax.transAxes)
    placed = []
    for tool, value in ends.items():
        pos = to_axes.transform((1024, value))[1]
        pos = max(pos, placed[-1] + 0.1) if placed else max(pos, 0.06)  # clear of the x axis
        placed.append(pos)
        ax.text(1024 * 1.04, pos, seconds_text(value), transform=at, va="center", fontsize=6, color=INK,
                weight="bold" if value == ends.min() else "normal")


def bars(axes, scan):
    """One panel per data type, each with its own y axis in seconds, the
    fastest bar's label in bold."""
    rows = largest[largest.scan == scan]
    for ax, dataset in zip(axes, DATASETS):
        data = rows[rows.dataset == dataset].set_index("tool").loc[list(TOOLS)]
        ax.bar(range(3), data.s, width=0.75, color=[c for _, c in TOOLS.values()], edgecolor="white", lw=0.8)
        for x, s in enumerate(data.s):
            ax.text(x, s, seconds_text(s), ha="center", va="bottom", fontsize=6, color=INK,
                    weight="bold" if s == data.s.min() else "normal")
        ax.margins(y=0.15)
        ax.set_xticks(range(3), [name for name, _ in TOOLS.values()])
        ax.set_title(f"{dataset}\n(n={data.n.iloc[0]:,})", loc="left")
    axes[0].set_ylabel("Wall time (s)")
    axes[0].annotate(SCANS[scan], (0, 1), xycoords="axes fraction", xytext=(0, 22), textcoords="offset points",
                     va="bottom", fontsize=7, weight="bold")


WIDTH, HEIGHT = 180 / 25.4, 5.9  # inches: the text width of an A4 page

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.titlesize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none",
})

# Rows (positions in inches): the legend, A across the page, the two B panels,
# then C (whole genome) and D (regions) side by side.
fig = plt.figure(figsize=(WIDTH, HEIGHT))
a = fig.add_axes([0.3 / WIDTH, 4.05 / HEIGHT, 6.6 / WIDTH, 1.45 / HEIGHT])
schematic(a, 6.6, 1.45)
row_b = fig.add_gridspec(1, 2, left=0.55 / WIDTH, right=6.6 / WIDTH, top=3.55 / HEIGHT, bottom=2.35 / HEIGHT,
                         wspace=0.3)
b = [fig.add_subplot(row_b[i]) for i in (0, 1)]
for ax, scan in zip(b, SCANS):
    scaling(ax, scan)
bottoms = {}
for scan, (left, right) in zip(SCANS, [(0.5, 3.35), (4.1, 6.95)]):
    grid = fig.add_gridspec(1, 3, left=left / WIDTH, right=right / WIDTH, top=1.45 / HEIGHT, bottom=0.3 / HEIGHT,
                            wspace=0.4)
    bottoms[scan] = [fig.add_subplot(grid[i]) for i in (0, 1, 2)]
    bars(bottoms[scan], scan)
letters = [letter(a, "A", dx=-4), letter(b[0], "B"),
           letter(bottoms["full_scan"][0], "C", dy=22), letter(bottoms["region_scan"][0], "D", dy=22)]
legend = fig.legend([Patch(color=c) for _, c in TOOLS.values()], [n for n, _ in TOOLS.values()],
                    loc="upper center", bbox_to_anchor=(0.5, 1), ncol=3, frameon=False)
# Cut the files tightly around what is drawn, so no empty canvas is left.
fig.savefig(ROOT / "figures" / "figure3.png", dpi=300, bbox_inches="tight", pad_inches=0.04)
fig.savefig(ROOT / "figures" / "figure3.svg", bbox_inches="tight", pad_inches=0.04)
PANELS.mkdir(exist_ok=True)
save_panels(fig, {"figure3-A": [a], "figure3-B": b, "figure3-C": bottoms["full_scan"],
                  "figure3-D": bottoms["region_scan"], "legend": [legend]}, letters)
