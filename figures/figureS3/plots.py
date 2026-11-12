"""Figure S3: wall time and peak memory of building the largest cohort of each data type."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results"
OUT = HERE / "plots"
plt.style.use(HERE.parent / "paper.mplstyle")
TOOLS = {"bigwig": ("BigWig", "#9CC9C0"), "d4": ("D4", "#A7BEDB"), "pbz": ("PBZ", "#C0304F")}
INK = "#444"


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
OUT.mkdir(exist_ok=True)
for ext in ("svg", "png"):
    fig.savefig(OUT / f"figureS3.{ext}")
