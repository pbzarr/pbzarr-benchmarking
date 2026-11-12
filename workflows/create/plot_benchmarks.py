#!/usr/bin/env python3
"""Plot aggregated benchmark and size TSVs with Seaborn."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, NullLocator
import pandas as pd
import seaborn as sns


LABELS = {"pbz": "PBZ", "d4": "D4", "bigwig": "BigWig"}
COLORS = {"PBZ": "#0072B2", "D4": "#D55E00", "BigWig": "#009E73"}
MARKERS = {"PBZ": "o", "D4": "s", "BigWig": "^"}
# Small multiplicative x offsets so replicate points of different tools don't overlap.
DODGE = {"PBZ": 0.94, "D4": 1.0, "BigWig": 1.06}
DATASETS = {
    "1kg-sr": "1000 Genomes 30x short read",
    "1kg-ont": "1000 Genomes ONT",
    "1kg-rna-ont": "1000 Genomes ONT direct RNA",
}
GIB = FuncFormatter(lambda value, _: f"{value:g}")


def panel(ax, data, metric, label, log=False):
    """Median line across replicates, with every replicate drawn as a point."""
    for tool in COLORS:
        rows = data[data.Tool == tool]
        if rows.empty:
            continue
        replicated = rows.groupby("n_samples").size().max() > 1
        dodge = DODGE[tool] if replicated else 1.0
        median = rows.groupby("n_samples")[metric].median()
        ax.plot(median.index * dodge, median.values, color=COLORS[tool],
                lw=2, marker=MARKERS[tool], ms=8, mec="white", mew=1.2, zorder=3)
        if replicated:
            ax.scatter(rows.n_samples * dodge, rows[metric], color=COLORS[tool], marker=MARKERS[tool],
                       s=22, alpha=0.45, lw=0, zorder=2)
    ticks = sorted(data.n_samples.unique())
    ax.set_xscale("log", base=2)
    ax.set_xticks(ticks, [str(t) for t in ticks])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlim(ticks[0] / 1.3, ticks[-1] * 1.3)
    ax.set(xlabel="Samples", ylabel=label)
    ax.grid(False, axis="x")
    if log:
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(GIB)
    else:
        ax.set_ylim(bottom=0)
    sns.despine(ax=ax)


def save_figure(fig, path, title):
    fig.suptitle(title, x=0.065, y=0.99, ha="left", fontsize=15, weight="bold")
    handles = [
        Line2D([], [], color=color, marker=MARKERS[tool], ms=8, mec="white", lw=2, label=tool)
        for tool, color in COLORS.items()
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.98, 0.995),
               ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.95), w_pad=2.5, h_pad=2.2)
    for extension in ("png", "pdf"):
        output = path.with_suffix(f".{extension}")
        fig.savefig(output, dpi=220, bbox_inches="tight", facecolor="white")
        print(output)
    plt.close(fig)


def main():
    root = Path(__file__).resolve().parents[2] / "results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("benchmarks", nargs="?", type=Path, default=root / "benchmarks.tsv")
    parser.add_argument("--sizes", type=Path, help="Defaults to sizes.tsv beside benchmarks")
    parser.add_argument("--outdir", type=Path, help="Defaults to plots/ beside benchmarks")
    args = parser.parse_args()
    sizes_path = args.sizes or args.benchmarks.parent / "sizes.tsv"
    outdir = args.outdir or args.benchmarks.parent / "plots"
    outdir.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(args.benchmarks, sep="\t")
    data = data[data.rule.isin(["create", "full_scan", "region_scan"])].copy()
    data["Tool"] = data.tool.map(LABELS)
    data["minutes"] = data.s / 60
    data["rss_gib"] = data.max_rss / 1024
    data["read_gib"] = data.io_in / 1024
    sizes = pd.read_csv(sizes_path, sep="\t") if sizes_path.exists() else None
    if sizes is not None:
        sizes["Tool"] = sizes.tool.map(LABELS)
        sizes["size_gib"] = sizes.size_kib / 1024**2

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2,
                  rc={"grid.color": "#E5E5E5", "axes.titleweight": "bold",
                      "pdf.fonttype": 42, "ps.fonttype": 42})
    for (dataset, threads), group in data.groupby(["dataset", "threads"], sort=True):
        threads = int(threads)
        name = DATASETS.get(dataset, dataset)
        creation = group[group.rule == "create"]
        disk = None if sizes is None else sizes[(sizes.dataset == dataset) & (sizes.threads == threads)]
        has_disk = disk is not None and not disk.empty
        if not creation.empty:
            fig, axes = plt.subplots(1, 3 if has_disk else 2, figsize=(13.5 if has_disk else 9, 4.6))
            panel(axes[0], creation, "minutes", "Wall time (min)")
            panel(axes[1], creation, "rss_gib", "Peak RSS (GiB)")
            axes[0].set_title("Creation time")
            axes[1].set_title("Creation memory")
            if has_disk:
                panel(axes[2], disk, "size_gib", "Allocated disk (GiB)")
                axes[2].set_title("Disk footprint")
                d4 = disk[disk.tool == "d4"].groupby("n_samples").size_gib.median()
                if not d4.empty:
                    axes[2].plot(d4.index, 2 * d4.values, ":", color=COLORS["D4"], lw=2)
                    axes[2].set_ylim(0, 1.12 * max(disk.size_gib.max(), 2 * d4.max()))
            save_figure(fig, outdir / f"{dataset}_t{threads}_creation",
                        f"{name} | Creation, {threads} threads")

        scans = [("full_scan", "Full scan"), ("region_scan", "Region scan")]
        if group.rule.isin([rule for rule, _ in scans]).any():
            fig, axes = plt.subplots(4, 2, figsize=(10, 13.5), sharex=True)
            for column, (rule, title) in enumerate(scans):
                subset = group[group.rule == rule]
                if subset.empty:
                    axes[0, column].set_visible(False)
                    continue
                panel(axes[0, column], subset, "s", "Wall time (s, log)", log=True)
                panel(axes[1, column], subset, "cpu_time", "CPU time (s, log)", log=True)
                panel(axes[2, column], subset, "rss_gib", "Peak RSS (GiB)")
                panel(axes[3, column], subset, "read_gib", "Read from storage (GiB)")
                axes[0, column].set_title(title)
                for row in (0, 1, 2):
                    axes[row, column].set_xlabel("")
            save_figure(fig, outdir / f"{dataset}_t{threads}_analysis",
                        f"{name} | Analysis, {threads} threads")


if __name__ == "__main__":
    main()
