#!/usr/bin/env python3
"""Figure 2, creation.

A: how each format builds a cohort depth file from N CRAMs.
B: short-read sample scaling: size on disk, wall time and peak memory.
C: size on disk of the largest cohort of each data type, with the percent of
   the input CRAMs on each bar.

Reads the builds.tsv files copied from Redwood into results/ and writes
figures/figure2.png and .svg (the whole figure), plus each panel as an SVG for
Illustrator in figures/panels/ (figure2-A, figure2-B, figure2-C, legend). The
figure is drawn at print size for an A4 page (180 mm wide, 6-7 pt text), so it
and the panels are placed at 100%. Run from the repo root:
    pixi run python figures/figure2.py
"""

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
# BigWig, then D4, then PBZ: the order the tools came in. PBZ is the vibrant
# color; the other two are pastel, with darker shades for their names.
TOOLS = {"bigwig": ("BigWig", "#9CC9C0"), "d4": ("D4", "#A7BEDB"), "pbz": ("PBZ", "#C0304F")}
DARK = {"bigwig": "#3F8073", "d4": "#4A6E9C", "pbz": "#C0304F"}
INK = "#444"
MONO = "Courier New"
COLUMNS = ["n", "tool", "read_bytes", "wall_s", "max_rss_mib", "out_bytes"]


def read(*patterns):
    files = [f for pattern in patterns for f in sorted(RESULTS.glob(pattern))]
    return pd.concat([pd.read_csv(f, sep="\t", usecols=COLUMNS) for f in files])


# Short read: PBZ and D4/BigWig ran as separate jobs. N = 256 has 3 builds.
sr = read("cohort-pbz/n*/b*/builds.tsv", "cohort-d4bw/n*/b*/builds.tsv")
sr["size_tb"] = sr.out_bytes / 1e12
sr["percent"] = 100 * sr.out_bytes / sr.read_bytes
sr["minutes"] = sr.wall_s / 60
sr["memory_gb"] = sr.max_rss_mib * 2**20 / 1e9

largest = pd.concat([
    sr[sr.n == 1024].assign(dataset="Short read"),
    read("cohort-ont/n1019/builds.tsv").assign(dataset="Long read"),
    read("cohort-rna/n60/builds.tsv").assign(dataset="Long-read RNA"),
])
largest["percent"] = 100 * largest.out_bytes / largest.read_bytes
DATASETS = list(largest.dataset.unique())


def size_text(n_bytes):
    if n_bytes >= 1e12:
        return f"{n_bytes / 1e12:.1f} TB"
    if n_bytes >= 1e10:
        return f"{n_bytes / 1e9:.0f} GB"
    return f"{n_bytes / 1e9:.1f} GB"


def percent_text(p):
    return f"{p:.1f}%" if p < 10 else f"{p:.0f}%"


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


def box(ax, x, y, text, tool=False, w=1.1, missing=False):
    """A file (white) or a tool (grey, monospace) centered on x, y. A missing
    file is a dashed grey outline with grey italic text."""
    ax.add_patch(FancyBboxPatch((x - w / 2, y - 0.19), w, 0.38, boxstyle="round,pad=0,rounding_size=0.08",
                                fc="#eeeeee" if tool else "white", ec="#999" if missing else INK, lw=0.6,
                                ls="--" if missing else "-", zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=5.5 if tool else 6, family=MONO if tool else None,
            color="#777" if missing else "black", style="italic" if missing else "normal", zorder=4)


def arrow(ax, start, end):
    ax.annotate("", end, start, zorder=2, arrowprops=dict(
        arrowstyle="-|>", mutation_scale=5, shrinkA=0, shrinkB=0, lw=0.6, color=INK))


def join(ax, xs, top, cx, bottom):
    """Pipeline-style arrows from each x at top into one box at cx, bottom: down,
    across at a right angle, then one arrowhead down into the box."""
    mid = bottom + 0.12  # below the dashed row line
    for x in xs:
        ax.plot([x, x, cx], [top, mid, mid], color=INK, lw=0.6, solid_joinstyle="miter", zorder=2)
    arrow(ax, (cx, mid), (cx, bottom))


def schematic(ax, width, height):
    """One column per tool: the N CRAMs (s1 to sN) at the top, the per-sample
    depth files they become, and the cohort file. PBZ goes straight from the
    CRAMs to the cohort file; BigWig has no merge, so no cohort file."""
    w, h = width / UNIT, height / UNIT
    ax.set(xlim=(0, w), ylim=(0, h), aspect="equal")
    ax.axis("off")
    head, raw, step1, profile, step2, cohort = h - 0.3, h - 0.95, h - 1.7, h - 2.45, h - 3.2, h - 3.95
    for y, label in [(raw, "Input\nalignments"), (step1, "Compute\ndepth"), (profile, "Per-sample\nfile"),
                     (step2, "Merge"), (cohort, "Cohort\ndepth")]:
        ax.text(0.05, y, label, va="center", fontsize=6, color=INK)
    # Faint dashed lines between the rows.
    rows = [raw, step1, profile, step2, cohort]
    for upper, lower in zip(rows, rows[1:]):
        ax.plot([0, w], [(upper + lower) / 2] * 2, color="#bbb", lw=0.4, ls=(0, (3, 2)), zorder=1)
    col = (w - 1.15) / 3
    for i, (tool, (name, color)) in enumerate(TOOLS.items()):
        cx = 1.15 + col * (i + 0.5)
        left, right = cx - 0.7, cx + 0.7
        ax.add_patch(FancyBboxPatch((cx - col / 2 + 0.06, 0), col - 0.12, h, boxstyle="round,pad=0,rounding_size=0.12",
                                    fc=color, alpha=0.25, lw=0, zorder=0))
        ax.text(cx, head, name, ha="center", va="center", fontsize=7, weight="bold", color=DARK[tool])
        for x, s in [(left, "s1"), (right, "sN")]:
            box(ax, x, raw, f"{s}.cram")
        ax.text(cx, raw, "…", ha="center", va="center", fontsize=6)
        if tool == "pbz":
            # PBZ computes depth straight into the cohort file: no per-sample files, no merge.
            join(ax, (left, right), raw - 0.19, cx, step1 + 0.19)
            box(ax, cx, step1, "pbz import", tool=True, w=1.3)
            arrow(ax, (cx, step1 - 0.19), (cx, cohort + 0.19))
            box(ax, cx, cohort, "cohort.pbz", w=1.3)
            continue
        maker, ext = {"bigwig": ("megadepth", "bw"), "d4": ("mosdepth", "d4")}[tool]
        for x, s in [(left, "s1"), (right, "sN")]:
            arrow(ax, (x, raw - 0.19), (x, step1 + 0.19))
            box(ax, x, step1, maker, tool=True, w=1.28)
            arrow(ax, (x, step1 - 0.19), (x, profile + 0.19))
            box(ax, x, profile, f"{s}.{ext}")
        ax.text(cx, profile, "…", ha="center", va="center", fontsize=6)
        if tool == "d4":  # BigWig has no merge tool, so no cohort file
            join(ax, (left, right), profile - 0.19, cx, step2 + 0.19)
            box(ax, cx, step2, "d4tools merge", tool=True, w=1.6)
            arrow(ax, (cx, step2 - 0.19), (cx, cohort + 0.19))
            box(ax, cx, cohort, "cohort.d4", w=1.3)


def scaling(ax, y, label, text, xmax=1400):
    """One line per tool over cohort size (the median of the builds), the three
    builds at N = 256 as a vertical range, and each tool's value at N = 1,024
    written past the line's end, the best (lowest) in bold."""
    for tool, (_, color) in TOOLS.items():
        rows = sr[sr.tool == tool]
        median = rows.groupby("n")[y].median()
        ax.plot(median.index, median.values, color=color, lw=1.6 if tool == "pbz" else 1, marker="o", ms=2.5, zorder=3)
        builds = rows[rows.n == 256][y]
        ax.vlines(256, builds.min(), builds.max(), color=color, lw=2.5, zorder=2)
    ax.set_ylabel(label)
    ax.set_xticks([0, 256, 512, 768, 1024])
    ax.set(xlim=(0, xmax), ylim=(0, None))

    # Spread the end labels so they are at least a line apart.
    ends = sr[sr.n == 1024].groupby("tool")[y].median().sort_values()
    to_axes = ax.transData + ax.transAxes.inverted()
    at = blended_transform_factory(ax.transData, ax.transAxes)
    placed = []
    for tool, value in ends.items():
        pos = to_axes.transform((1024, value))[1]
        pos = max(pos, placed[-1] + 0.1) if placed else pos
        placed.append(pos)
        ax.text(1024 * 1.04, pos, text(tool, value), transform=at, va="center", fontsize=6, color=INK,
                weight="bold" if value == ends.min() else "normal")


def panel_b(axes):
    def size(tool, tb):
        percent = sr[(sr.n == 1024) & (sr.tool == tool)].percent.median()
        return f"{tb:.2g} TB ({percent:.0f}% of input)"
    scaling(axes[0], "size_tb", "Size on disk (TB)", size, xmax=1750)
    scaling(axes[1], "minutes", "Wall time (min)", lambda t, v: f"{v:.0f} min")
    scaling(axes[2], "memory_gb", "Peak memory (GB)", lambda t, v: f"{v:,.0f} GB")
    # Memory levels off at N = 256: the node has 256 cores, so at most 256 jobs run at once.
    axes[2].axvline(256, color="#999", lw=0.6, ls="--", zorder=1)
    axes[2].set_xlabel("Samples")


def panel_c(axes):
    """One panel per data type, each with its own y axis, labeled with the size
    and the percent of the input CRAMs."""
    for ax, dataset in zip(axes, DATASETS):
        rows = largest[largest.dataset == dataset].set_index("tool").loc[list(TOOLS)]
        scale, unit = (1e12, "TB") if rows.out_bytes.max() >= 1e12 else (1e9, "GB")
        ax.bar(range(3), rows.out_bytes / scale, width=0.75, color=[c for _, c in TOOLS.values()],
               edgecolor="white", lw=0.8)
        for x, (tool, row) in enumerate(rows.iterrows()):
            ax.text(x, row.out_bytes / scale, f"{size_text(row.out_bytes)}\n{percent_text(row.percent)}",
                    ha="center", va="bottom", fontsize=6, color=INK, linespacing=1,
                    weight="bold" if tool == "pbz" else "normal")
        ax.margins(y=0.3)
        ax.set_ylabel(f"Size on disk ({unit})")
        ax.set_xticks(range(3), [name for name, _ in TOOLS.values()])
        ax.set_title(f"{dataset} (n={rows.n.iloc[0]:,})\n{size_text(rows.read_bytes.iloc[0])} input", loc="left")


WIDTH, HEIGHT = 180 / 25.4, 5.5  # inches: the text width of an A4 page

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.titlesize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none",
})

# Two columns (positions in inches): A above C on the left, the three B panels
# stacked on the right, and the legend in one row on top.
fig = plt.figure(figsize=(WIDTH, HEIGHT))
a = fig.add_axes([0.3 / WIDTH, 3.0 / HEIGHT, 4.0 / WIDTH, 2.05 / HEIGHT])
schematic(a, 4.0, 2.05)
column_b = fig.add_gridspec(3, 1, left=4.95 / WIDTH, right=6.95 / WIDTH, top=4.95 / HEIGHT,
                            bottom=0.45 / HEIGHT, hspace=0.25)
b = [fig.add_subplot(column_b[i]) for i in (0, 1, 2)]
panel_b(b)
fig.align_ylabels(b)
grid = fig.add_gridspec(1, 3, left=0.45 / WIDTH, right=4.25 / WIDTH, top=2.3 / HEIGHT, bottom=0.3 / HEIGHT,
                        wspace=0.6)
c = [fig.add_subplot(grid[i]) for i in (0, 1, 2)]
panel_c(c)
letters = [letter(a, "A", dx=-4), letter(b[0], "B", dx=-30), letter(c[0], "C", dy=24)]
legend = fig.legend([Patch(color=c) for _, c in TOOLS.values()], [n for n, _ in TOOLS.values()],
                    loc="upper center", bbox_to_anchor=(0.5, 1), ncol=3, frameon=False)
# Cut the files tightly around what is drawn, so no empty canvas is left.
fig.savefig(ROOT / "figures" / "figure2.png", dpi=300, bbox_inches="tight", pad_inches=0.04)
fig.savefig(ROOT / "figures" / "figure2.svg", bbox_inches="tight", pad_inches=0.04)
PANELS.mkdir(exist_ok=True)
save_panels(fig, {"figure2-A": [a], "figure2-B": b, "figure2-C": c, "legend": [legend]}, letters)
