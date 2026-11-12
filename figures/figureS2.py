#!/usr/bin/env python3
"""Figure S2, what a short-read build costs.

A: total CPU time of the build over cohort size.
B: CPU time per sample over cohort size.
C: D4's wall time, split into the per-sample mosdepth phase and d4tools merge.

Lines and bars are the median of the builds at each N (N = 256 has 3 builds).
Reads the builds.tsv and merge.tsv files in results/ and writes
figures/figureS2.png and .svg at print size (180 mm wide). Run from the repo root:
    pixi run python figures/figureS2.py
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
D4_DARK = "#4A6E9C"
INK = "#444"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.titlesize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none",
})


def read(pattern):
    return pd.concat([pd.read_csv(f, sep="\t") for f in sorted(RESULTS.glob(pattern))])


# Short read: PBZ and D4/BigWig ran as separate jobs.
builds = pd.concat([read("cohort-pbz/n*/b*/builds.tsv"), read("cohort-d4bw/n*/b*/builds.tsv")])
builds["cpu_hours"] = builds.cpu_s / 3600
builds["cpu_minutes_per_sample"] = builds.cpu_s / 60 / builds.n

# D4: the merge is timed on its own; the rest of the build is the per-sample phase.
merges = read("cohort-d4bw/n*/b*/merge.tsv").rename(columns={"wall_s": "merge_s"})
d4 = builds[builds.tool == "d4"].merge(merges[["n", "build", "merge_s"]], on=["n", "build"])
d4 = d4.groupby("n")[["wall_s", "merge_s"]].median() / 60  # minutes
d4["per_sample"] = d4.wall_s - d4.merge_s


def lines(ax, column, label):
    """One line per tool over cohort size."""
    for tool, (_, color) in TOOLS.items():
        median = builds[builds.tool == tool].groupby("n")[column].median()
        ax.plot(median.index, median.values, color=color, lw=1.6 if tool == "pbz" else 1, marker="o", ms=2.5)
    ax.set(xlabel="Samples", ylabel=label, ylim=(0, None), xticks=[0, 256, 512, 768, 1024])


fig, (a, b, c) = plt.subplots(1, 3, figsize=(180 / 25.4, 2.4), layout="constrained")

lines(a, "cpu_hours", "CPU time (core-hours)")
a.set_title("Total CPU time", loc="left", weight="bold")

lines(b, "cpu_minutes_per_sample", "CPU time per sample (core-minutes)")
b.set_title("CPU time per sample", loc="left", weight="bold")

x = range(len(d4))
c.bar(x, d4.per_sample, width=0.75, color=TOOLS["d4"][1], edgecolor="white", lw=0.8)
c.bar(x, d4.merge_s, width=0.75, bottom=d4.per_sample, color=D4_DARK, edgecolor="white", lw=0.8)
for i, row in enumerate(d4.itertuples()):
    c.text(i, row.wall_s, f"{100 * row.merge_s / row.wall_s:.0f}%", ha="center", va="bottom", fontsize=6, color=INK)
c.set_xticks(x, d4.index)
c.set(xlabel="Samples", ylabel="Wall time (min)")
c.margins(y=0.12)
c.set_title("D4 build (label: merge share)", loc="left", weight="bold")
c.legend([Patch(color=TOOLS["d4"][1]), Patch(color=D4_DARK)], ["mosdepth, per sample", "d4tools merge"],
         frameon=False, loc="upper left", fontsize=6)

for ax, text in zip((a, b, c), "ABC"):
    ax.text(-0.2, 1.08, text, transform=ax.transAxes, fontsize=10, weight="bold")
fig.legend([Patch(color=color) for _, color in TOOLS.values()], [name for name, _ in TOOLS.values()],
           loc="outside upper center", ncol=3, frameon=False)
fig.savefig(ROOT / "figures" / "figureS2.png", dpi=300)
fig.savefig(ROOT / "figures" / "figureS2.svg")
