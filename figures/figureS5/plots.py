"""Figure S5: CPU time and peak memory of the short-read scans (log scale)."""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results"
OUT = HERE / "plots"
plt.style.use(HERE.parent / "paper.mplstyle")
TOOLS = {"bigwig": ("BigWig", "#9CC9C0"), "d4": ("D4", "#A7BEDB"), "pbz": ("PBZ", "#C0304F")}
SCANS = {"full_scan": "Mean depth, whole genome", "region_scan": "Mean depth, 10,000 random 10 kb regions"}
MEASURES = {"cpu_hours": "CPU time (core-hours)", "memory_gb": "Peak memory (GB)"}


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
OUT.mkdir(exist_ok=True)
for ext in ("svg", "png"):
    fig.savefig(OUT / f"figureS5.{ext}")
