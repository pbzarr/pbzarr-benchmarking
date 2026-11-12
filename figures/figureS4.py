#!/usr/bin/env python3
"""Figure S4, how much data each scan reads from disk (short read).

A: data read by the whole-genome scan over cohort size.
B: data read by the region scan over cohort size (log scale).
C: the share of the cohort's files that each scan reads at N = 1,024.

Every scan starts with a cold page cache, so the data read is what the scan
needs from disk. Values are the median of 3 trials.

Reads the scans/*.tsv and builds.tsv files in results/ and writes
figures/figureS4.png and .svg at print size (180 mm wide). Run from the repo root:
    pixi run python figures/figureS4.py
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
SCANS = {"full_scan": "Whole genome", "region_scan": "10,000 regions"}
INK = "#444"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.titlesize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none",
})

# One row per scan: N, tool, scan and the GB read (io_in is in MiB).
rows = []
for folder in ("cohort-pbz", "cohort-d4bw"):
    for f in sorted(RESULTS.glob(f"{folder}/n*/b1/scans/*_scan.tsv")):
        tool, scan = f.stem.split(".")
        n = int(re.search(r"/n(\d+)/", str(f)).group(1))
        read_gb = pd.read_csv(f, sep="\t").io_in.median() * 2**20 / 1e9
        rows.append({"n": n, "tool": tool, "scan": scan, "read_gb": read_gb})
scans = pd.DataFrame(rows)

# The size of each cohort's files, to get the share that a scan reads.
builds = pd.concat(pd.read_csv(f, sep="\t") for folder in ("cohort-pbz", "cohort-d4bw")
                   for f in sorted(RESULTS.glob(f"{folder}/n*/b1/builds.tsv")))
scans = scans.merge(builds[["n", "tool", "out_bytes"]], on=["n", "tool"])
scans["percent"] = 100 * scans.read_gb / (scans.out_bytes / 1e9)


def lines(ax, scan, label, scale=1):
    """One line per tool over cohort size."""
    for tool, (_, color) in TOOLS.items():
        line = scans[(scans.scan == scan) & (scans.tool == tool)].sort_values("n")
        ax.plot(line.n, line.read_gb / scale, color=color, lw=1.6 if tool == "pbz" else 1, marker="o", ms=2.5)
    ax.set(xlabel="Samples", ylabel=label, xticks=[0, 256, 512, 768, 1024])
    ax.set_title(SCANS[scan], loc="left", weight="bold")


fig, (a, b, c) = plt.subplots(1, 3, figsize=(180 / 25.4, 2.4), layout="constrained")

lines(a, "full_scan", "Data read (TB)", scale=1000)
a.set_ylim(0, None)

lines(b, "region_scan", "Data read (GB)")
b.set_yscale("log")

# Two groups of three bars: the share of the files read by each scan.
largest = scans[scans.n == 1024]
for group, scan in enumerate(SCANS):
    for i, (tool, (_, color)) in enumerate(TOOLS.items()):
        percent = largest[(largest.scan == scan) & (largest.tool == tool)].percent.iloc[0]
        x = group * 4 + i
        c.bar(x, percent, width=0.9, color=color, edgecolor="white", lw=0.8)
        c.text(x, percent, f"{percent:.0f}%", ha="center", va="bottom", fontsize=6, color=INK,
               weight="bold" if tool == "pbz" else "normal")
c.set_xticks([1, 5], SCANS.values())
c.set(ylabel="Share of the files read (%)", ylim=(0, 115))
c.set_title("N = 1,024", loc="left", weight="bold")

for ax, text in zip((a, b, c), "ABC"):
    ax.text(-0.2, 1.08, text, transform=ax.transAxes, fontsize=10, weight="bold")
fig.legend([Patch(color=color) for _, color in TOOLS.values()], [name for name, _ in TOOLS.values()],
           loc="outside upper center", ncol=3, frameon=False)
fig.savefig(ROOT / "figures" / "figureS4.png", dpi=300)
fig.savefig(ROOT / "figures" / "figureS4.svg")
