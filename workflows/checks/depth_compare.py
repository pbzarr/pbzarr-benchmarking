"""Compare per-base depth on one contig across tool outputs.

Reads every <label>.d4, <label>.bw and <label>.pbz in OUTDIR, loads CONTIG from each,
and prints how many bases differ between each pair. For each label it also writes
<label>.vs.<baseline>.bed: the runs of bases where it differs from the baseline, with the
largest difference in the run, so the reads behind them can be looked at.

usage: python depth_compare.py OUTDIR CONTIG BASELINE
"""

import io
import subprocess
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import pbzarr as pbz
import pyBigWig

outdir, contig, baseline = Path(sys.argv[1]), sys.argv[2], sys.argv[3]


def load(path):
    if path.suffix == ".d4":
        # d4tools view, not pyd4: our pyd4 build misreads files that d4tools create writes.
        # view prints every interval of the contig, zeros included.
        view = subprocess.run(["d4tools", "view", str(path), contig], capture_output=True, check=True)
        bg = pd.read_csv(io.BytesIO(view.stdout), sep="\t", header=None, usecols=[1, 2, 3]).to_numpy()
        return np.repeat(bg[:, 2], bg[:, 1] - bg[:, 0]).astype(np.int64)
    if path.suffix == ".bw":
        with pyBigWig.open(str(path)) as bw:
            values = bw.values(contig, 0, bw.chroms()[contig], numpy=True)
        return np.nan_to_num(values).astype(np.int64)  # no record = 0
    ds = pbz.open(str(path))["depth"].to_dataset(inherit=False)
    return ds.pbz.region(contig).values.astype(np.int64)


depth = {p.stem: load(p) for p in sorted(outdir.iterdir()) if p.suffix in {".d4", ".bw", ".pbz"}}
base = depth[baseline]
covered = base > 0

print(f"{contig}: {len(base):,} bases, {covered.sum():,} covered in {baseline}\n")
print(f"{'tool':<22}{'mean depth':>12}{'max':>8}{'bases > 0':>14}")
for label, d in depth.items():
    print(f"{label:<22}{d.mean():>12.3f}{d.max():>8}{(d > 0).sum():>14,}")

print(f"\n{'pair':<44}{'bases differ':>14}{'% of covered':>14}{'a > b':>12}{'a < b':>12}{'max |a-b|':>11}")
for a, b in combinations(depth, 2):
    diff = depth[a] - depth[b]
    n = np.count_nonzero(diff)
    print(f"{a + ' vs ' + b:<44}{n:>14,}{100 * n / covered.sum():>13.3f}%"
          f"{(diff > 0).sum():>12,}{(diff < 0).sum():>12,}{np.abs(diff).max():>11}")

for label, d in depth.items():
    if label == baseline:
        continue
    diff = d - base
    differs = np.flatnonzero(diff)
    if not len(differs):
        continue
    # runs of consecutive differing bases -> BED intervals
    breaks = np.flatnonzero(np.diff(differs) > 1) + 1
    with open(outdir / f"{label}.vs.{baseline}.bed", "w") as bed:
        for run in np.split(differs, breaks):
            biggest = run[np.argmax(np.abs(diff[run]))]
            bed.write(f"{contig}\t{run[0]}\t{run[-1] + 1}\t{diff[biggest]}\t{base[biggest]}\t{d[biggest]}\n")
