"""Figure 3 plots: B, short-read scan scaling; C and D, whole-genome and region scans
of the largest cohort per data type; the legend. Each scan is the median of 3 cold-cache trials."""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results"
OUT = HERE / "plots"
plt.style.use(HERE.parent / "paper.mplstyle")

TOOLS = {"bigwig": ("BigWig", "#9CC9C0"), "d4": ("D4", "#A7BEDB"), "pbz": ("PBZ", "#C0304F")}
INK = "#444"
SCANS = {"full_scan": "Calculating mean coverage of the whole genome",
         "region_scan": "Calculating mean coverage of 10,000 random regions"}


def read(pattern):
    rows = []
    for f in sorted(RESULTS.glob(pattern)):
        tool, scan = f.stem.split(".")
        n = int(re.search(r"/n(\d+)/", str(f)).group(1))
        rows.append({"n": n, "tool": tool, "scan": scan, "s": pd.read_csv(f, sep="\t").s.median()})
    return pd.DataFrame(rows)


sr = pd.concat([read("cohort-pbz/n*/b1/scans/*_scan.tsv"), read("cohort-d4bw/n*/b1/scans/*_scan.tsv")])

largest = pd.concat([
    sr[sr.n == 1024].assign(dataset="Short read"),
    read("cohort-ont/n1019/scans/*_scan.tsv").assign(dataset="Long read"),
    read("cohort-rna/n60/scans/*_scan.tsv").assign(dataset="Long-read RNA"),
])


def seconds_text(s):
    return f"{s:,.0f} s" if s >= 10 else f"{s:.1f} s"


def save(fig, name):
    OUT.mkdir(exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"{name}.{ext}")


def scaling(ax, scan):
    rows = sr[sr.scan == scan]
    for tool, (_, color) in TOOLS.items():
        line = rows[rows.tool == tool].sort_values("n")
        ax.plot(line.n, line.s, color=color, lw=1.6 if tool == "pbz" else 1, marker="o", ms=2.5, zorder=3)
    ax.set_title(SCANS[scan], loc="left", weight="bold")
    ax.set(xlabel="Samples", ylabel="Wall time (s)", xlim=(0, 1400), ylim=(0, None))
    ax.set_xticks([0, 256, 512, 768, 1024])

    # Time at N = 1,024 past each line's end, at least a line apart, the fastest in bold.
    ends = rows[rows.n == 1024].set_index("tool").s.sort_values()
    to_axes = ax.transData + ax.transAxes.inverted()
    at = blended_transform_factory(ax.transData, ax.transAxes)
    placed = []
    for value in ends:
        pos = to_axes.transform((1024, value))[1]
        pos = max(pos, placed[-1] + 0.1) if placed else max(pos, 0.06)
        placed.append(pos)
        ax.text(1024 * 1.04, pos, seconds_text(value), transform=at, va="center", fontsize=6, color=INK,
                weight="bold" if value == ends.min() else "normal")


def bars(scan, name):
    fig, axes = plt.subplots(1, 3, figsize=(3.34, 1.81), layout="constrained")
    rows = largest[largest.scan == scan]
    for ax, (dataset, data) in zip(axes, rows.groupby("dataset", sort=False)):
        data = data.set_index("tool").loc[list(TOOLS)]
        ax.bar(range(3), data.s, width=0.75, color=[c for _, c in TOOLS.values()], edgecolor="white", lw=0.8)
        for x, s in enumerate(data.s):
            ax.text(x, s, seconds_text(s), ha="center", va="bottom", fontsize=6, color=INK,
                    weight="bold" if s == data.s.min() else "normal")
        ax.margins(y=0.15)
        ax.set_xticks(range(3), [name for name, _ in TOOLS.values()])
        ax.set_title(f"{dataset}\n(n={data.n.iloc[0]:,})", loc="left")
    axes[0].set_ylabel("Wall time (s)")
    fig.suptitle(SCANS[scan], x=0.02, ha="left", fontsize=7, weight="bold")
    save(fig, name)


fig, axes = plt.subplots(1, 2, figsize=(6.54, 1.76), layout="constrained")
for ax, scan in zip(axes, SCANS):
    scaling(ax, scan)
save(fig, "B")

bars("full_scan", "C")
bars("region_scan", "D")

fig = plt.figure(figsize=(2.03, 0.27))
fig.legend([Patch(color=c) for _, c in TOOLS.values()], [n for n, _ in TOOLS.values()],
           loc="center", ncol=3, frameon=False)
save(fig, "legend")
