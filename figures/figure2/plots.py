"""Figure 2 plots: B, short-read scaling; C, size of the largest cohort per data type; the legend."""

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
COLUMNS = ["n", "tool", "read_bytes", "wall_s", "max_rss_mib", "out_bytes"]


def read(*patterns):
    files = [f for pattern in patterns for f in sorted(RESULTS.glob(pattern))]
    return pd.concat([pd.read_csv(f, sep="\t", usecols=COLUMNS) for f in files])


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


def size_text(n_bytes):
    if n_bytes >= 1e12:
        return f"{n_bytes / 1e12:.1f} TB"
    if n_bytes >= 1e10:
        return f"{n_bytes / 1e9:.0f} GB"
    return f"{n_bytes / 1e9:.1f} GB"


def percent_text(p):
    return f"{p:.1f}%" if p < 10 else f"{p:.0f}%"


def save(fig, name):
    OUT.mkdir(exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"{name}.{ext}")


def scaling(ax, y, label, text, xmax=1400):
    for tool, (_, color) in TOOLS.items():
        rows = sr[sr.tool == tool]
        median = rows.groupby("n")[y].median()
        ax.plot(median.index, median.values, color=color, lw=1.6 if tool == "pbz" else 1, marker="o", ms=2.5, zorder=3)
        repeats = rows[rows.n == 256][y]
        ax.vlines(256, repeats.min(), repeats.max(), color=color, lw=2.5, zorder=2)
    ax.set_ylabel(label)
    ax.set_xticks([0, 256, 512, 768, 1024])
    ax.set(xlim=(0, xmax), ylim=(0, None))

    # Value at N = 1,024 past each line's end, at least a line apart, the lowest in bold.
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


def size_label(tool, tb):
    percent = sr[(sr.n == 1024) & (sr.tool == tool)].percent.median()
    return f"{tb:.2g} TB ({percent:.0f}% of input)"


fig, (size, time, memory) = plt.subplots(3, 1, figsize=(2.49, 4.91), layout="constrained")
scaling(size, "size_tb", "Size on disk (TB)", size_label, xmax=1750)
scaling(time, "minutes", "Wall time (min)", lambda t, v: f"{v:.0f} min")
scaling(memory, "memory_gb", "Peak memory (GB)", lambda t, v: f"{v:,.0f} GB")
memory.axvline(256, color="#999", lw=0.6, ls="--", zorder=1)
memory.set_xlabel("Samples")
fig.align_ylabels()
save(fig, "B")

fig, axes = plt.subplots(1, 3, figsize=(4.31, 2.53), layout="constrained")
for ax, (dataset, rows) in zip(axes, largest.groupby("dataset", sort=False)):
    rows = rows.set_index("tool").loc[list(TOOLS)]
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
save(fig, "C")

fig = plt.figure(figsize=(2.03, 0.27))
fig.legend([Patch(color=c) for _, c in TOOLS.values()], [n for n, _ in TOOLS.values()],
           loc="center", ncol=3, frameon=False)
save(fig, "legend")
