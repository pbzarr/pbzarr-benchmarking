"""Figure S2: D4's short-read build time, split into mosdepth (per sample) and d4tools merge."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1] / "results"
OUT = HERE / "plots"
plt.style.use(HERE.parent / "paper.mplstyle")
D4 = "#A7BEDB"
D4_DARK = "#4A6E9C"
INK = "#444"


def read(pattern):
    return pd.concat([pd.read_csv(f, sep="\t") for f in sorted(RESULTS.glob(pattern))])


# The merge is timed on its own; the rest of the build is the per-sample phase.
builds = read("cohort-d4bw/n*/b*/builds.tsv")
merges = read("cohort-d4bw/n*/b*/merge.tsv").rename(columns={"wall_s": "merge_s"})
d4 = builds[builds.tool == "d4"].merge(merges[["n", "build", "merge_s"]], on=["n", "build"])
d4 = d4.groupby("n")[["wall_s", "merge_s"]].median() / 60  # minutes
d4["per_sample"] = d4.wall_s - d4.merge_s

fig, ax = plt.subplots(figsize=(90 / 25.4, 2.2), layout="constrained")
x = range(len(d4))
ax.bar(x, d4.per_sample, width=0.75, color=D4, edgecolor="white", lw=0.4)
ax.bar(x, d4.merge_s, width=0.75, bottom=d4.per_sample, color=D4_DARK, edgecolor="white", lw=0.4)
for i, row in enumerate(d4.itertuples()):
    ax.text(i, row.wall_s, f"{100 * row.merge_s / row.wall_s:.0f}%", ha="center", va="bottom", fontsize=6, color=INK)
ax.set_xticks(x, [f"{n:,}" for n in d4.index])
ax.set(xlabel="Samples", ylabel="Wall time (min)")
ax.margins(y=0.12)
ax.set_title("D4 build time", loc="left", weight="bold")
ax.legend([Patch(color=D4), Patch(color=D4_DARK)], ["mosdepth, per sample", "d4tools merge"],
          frameon=False, loc="upper left", fontsize=6)

OUT.mkdir(exist_ok=True)
for ext in ("svg", "png"):
    fig.savefig(OUT / f"figureS2.{ext}")
# Illustrator opens the PDF-format .ai with text kept as editable Arial.
fig.savefig(OUT / "figureS2.ai", format="pdf")
