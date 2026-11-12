#!/usr/bin/env python3
"""Figure S5, CPU time and peak memory of the scans (short read).

Top row: CPU time over cohort size. Bottom row: peak memory over cohort size.
Left: the whole-genome scan. Right: the region scan. All y axes are log scale.

Values are the median of 3 cold-cache trials. Peak memory is the resident
memory summed over all processes of the scan; BigWig runs one process per file,
up to 256 at once.

Reads the scans/*.tsv files in results/ and writes figures/figureS5.png and
.svg at print size (180 mm wide). Run from the repo root:
    pixi run python figures/figureS5.py
"""

import re
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
SCANS = {"full_scan": "Mean depth, whole genome", "region_scan": "Mean depth, 10,000 random 10 kb regions"}
MEASURES = {"cpu_hours": "CPU time (core-hours)", "memory_gb": "Peak memory (GB)"}

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.titlesize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none",
})

# One row per scan: N, tool, scan, and the median CPU time and peak memory
# (cpu_time is in seconds, max_rss in MiB).
rows = []
for folder in ("cohort-pbz", "cohort-d4bw"):
    for f in sorted(RESULTS.glob(f"{folder}/n*/b1/scans/*_scan.tsv")):
        tool, scan = f.stem.split(".")
        n = int(re.search(r"/n(\d+)/", str(f)).group(1))
        trials = pd.read_csv(f, sep="\t")
        rows.append({"n": n, "tool": tool, "scan": scan,
                     "cpu_hours": trials.cpu_time.median() / 3600,
                     "memory_gb": trials.max_rss.median() * 2**20 / 1e9})
scans = pd.DataFrame(rows)

fig, axes = plt.subplots(2, 2, figsize=(180 / 25.4, 4.0), layout="constrained", sharex=True)
for row, (measure, label) in enumerate(MEASURES.items()):
    for column, scan in enumerate(SCANS):
        ax = axes[row, column]
        for tool, (_, color) in TOOLS.items():
            line = scans[(scans.scan == scan) & (scans.tool == tool)].sort_values("n")
            ax.plot(line.n, line[measure], color=color, lw=1.6 if tool == "pbz" else 1, marker="o", ms=2.5)
        ax.set_yscale("log")
        ax.set_ylabel(label)
        ax.set_xticks([0, 256, 512, 768, 1024])
for ax, title in zip(axes[0], SCANS.values()):
    ax.set_title(title, loc="left", weight="bold")
for ax in axes[1]:
    ax.set_xlabel("Samples")

for ax, text in zip(axes[:, 0], "AB"):
    ax.text(-0.15, 1.06, text, transform=ax.transAxes, fontsize=10, weight="bold")
fig.legend([Patch(color=color) for _, color in TOOLS.values()], [name for name, _ in TOOLS.values()],
           loc="outside upper center", ncol=3, frameon=False)
fig.savefig(ROOT / "figures" / "figureS5.png", dpi=300)
fig.savefig(ROOT / "figures" / "figureS5.svg")
