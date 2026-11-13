"""Figure S3: single-sample build time, peak memory and scan times by thread count, for the median-size sample of each data type."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results"
OUT = HERE / "plots"
plt.style.use(HERE.parent / "paper.mplstyle")

TOOLS = {
    "megadepth": ("megadepth (BigWig)", "#9CC9C0"),
    "mosdepth": ("mosdepth (D4)", "#A7BEDB"),
    "pbz": ("pbz import (PBZ)", "#C0304F"),
}
DATASETS = {"single": "Short read", "single-ont": "Long read", "single-rna": "Long-read RNA"}
THREADS = [1, 4, 16, 64]


def runs(folder):
    """Every build run: tool, threads, wall_s, memory_gb."""
    files = sorted((RESULTS / folder / "runs").glob("*/*/t*/r*.tsv"))
    data = pd.concat([pd.read_csv(f, sep="\t") for f in files])
    data["memory_gb"] = data.max_rss_kb * 1024 / 1e9
    return data


def scans(folder):
    """Every scan trial, from scans/{tool}/{sample}/t{threads}/{scan}.tsv: tool, threads, scan, s."""
    trials = []
    for f in sorted((RESULTS / folder / "scans").glob("*/*/t*/*.tsv")):
        tool, _, threads, scan = f.with_suffix("").parts[-4:]
        trials.append(pd.read_csv(f, sep="\t").assign(tool=tool, threads=int(threads[1:]), scan=scan))
    return pd.concat(trials)


def dots(ax, data, column, label):
    """The median of the repeats as a dot joined across thread counts, their minimum to maximum as a bar.
    Tools sit side by side."""
    for offset, (tool, (_, color)) in zip([0.8, 1, 1.25], TOOLS.items()):
        stats = data[data.tool == tool].groupby("threads")[column].agg(["median", "min", "max"])
        median, low, high = (stats[c].to_numpy() for c in ("median", "min", "max"))
        ax.errorbar(stats.index.to_numpy() * offset, median, yerr=[median - low, high - median],
                    fmt="o-", lw=1.6 if tool == "pbz" else 1, ms=2.5, color=color, elinewidth=0.6, capsize=1.5, capthick=0.6)
    ax.set_xscale("log", base=2)
    ax.set_xticks(THREADS, THREADS)
    ax.minorticks_off()
    ax.set(xlabel="Threads", ylabel=label, ylim=(0, None))


fig, axes = plt.subplots(4, 3, figsize=(180 / 25.4, 7.4), layout="constrained")
for column, (folder, dataset) in enumerate(DATASETS.items()):
    builds = runs(folder)
    dots(axes[0, column], builds, "wall_s", "Wall time (s)")
    dots(axes[1, column], builds, "memory_gb", "Peak memory (GB)")
    trials = scans(folder)
    dots(axes[2, column], trials[trials.scan == "full_scan"], "s", "Full scan (s)")
    dots(axes[3, column], trials[trials.scan == "region_scan"], "s", "Region scan (s)")
    axes[0, column].set_title(dataset, loc="left", weight="bold")

for ax, text in zip(axes[:, 0], "ABCD"):
    ax.text(-0.25, 1.08, text, transform=ax.transAxes, fontsize=10, weight="bold")
fig.legend([Line2D([], [], color=color, marker="o", ms=2.5) for _, color in TOOLS.values()],
           [name for name, _ in TOOLS.values()], loc="outside upper center", ncol=3, frameon=False)
OUT.mkdir(exist_ok=True)
for ext in ("svg", "png"):
    fig.savefig(OUT / f"figureS3.{ext}")
# Illustrator opens the PDF-format .ai with text kept as editable Arial.
fig.savefig(OUT / "figureS3.ai", format="pdf")
