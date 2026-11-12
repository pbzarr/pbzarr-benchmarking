#!/usr/bin/env python3
"""Figure S3, building the largest cohort of each data type.

Top row: wall time of the build. Bottom row: peak memory. One column per data
type, each with its own y axis. The lowest bar's label is bold.

Peak memory is the resident memory summed over all processes of the build. For
D4 it includes the output file that each mosdepth process maps into memory.

Reads the builds.tsv files in results/ and writes figures/figureS3.png and .svg
at print size (180 mm wide). Run from the repo root:
    pixi run python figures/figureS3.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
# As in figures 2 and 3: BigWig, D4, PBZ, with PBZ the vibrant color.
TOOLS = {"bigwig": ("BigWig", "#9CC9C0"), "d4": ("D4", "#A7BEDB"), "pbz": ("PBZ", "#C0304F")}
INK = "#444"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.titlesize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none",
})


def read(*files):
    return pd.concat([pd.read_csv(RESULTS / f, sep="\t") for f in files])


# The same builds as figure 2C.
DATASETS = {
    "Short read": read("cohort-pbz/n1024/b1/builds.tsv", "cohort-d4bw/n1024/b1/builds.tsv"),
    "Long read": read("cohort-ont/n1019/builds.tsv"),
    "Long-read RNA": read("cohort-rna/n60/builds.tsv"),
}


def bars(ax, values, label):
    """One bar per tool, with its value on top."""
    ax.bar(range(3), values, width=0.75, color=[color for _, color in TOOLS.values()], edgecolor="white", lw=0.8)
    for x, value in enumerate(values):
        ax.text(x, value, f"{value:,.0f}" if value >= 10 else f"{value:.1f}", ha="center", va="bottom",
                fontsize=6, color=INK, weight="bold" if value == values.min() else "normal")
    ax.set_xticks(range(3), [name for name, _ in TOOLS.values()])
    ax.set_ylabel(label)
    ax.margins(y=0.15)


fig, axes = plt.subplots(2, 3, figsize=(180 / 25.4, 3.8), layout="constrained")
for column, (dataset, builds) in enumerate(DATASETS.items()):
    builds = builds.set_index("tool").loc[list(TOOLS)]
    bars(axes[0, column], builds.wall_s / 60, "Wall time (min)")
    bars(axes[1, column], builds.max_rss_mib * 2**20 / 1e9, "Peak memory (GB)")
    axes[0, column].set_title(f"{dataset} (n={builds.n.iloc[0]:,})", loc="left", weight="bold")

for ax, text in zip(axes[:, 0], "AB"):
    ax.text(-0.25, 1.08, text, transform=ax.transAxes, fontsize=10, weight="bold")
fig.legend([Patch(color=color) for _, color in TOOLS.values()], [name for name, _ in TOOLS.values()],
           loc="outside upper center", ncol=3, frameon=False)
fig.savefig(ROOT / "figures" / "figureS3.png", dpi=300)
fig.savefig(ROOT / "figures" / "figureS3.svg")
