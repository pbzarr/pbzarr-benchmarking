"""Where do the tools' per-base depths agree? Compare one contig, base by base,
for every pair of tool outputs, and write one row per pair.

A base counts as covered if either tool of the pair gives it a depth above 0.
differ_pct is the share of those bases where the two depths are not equal.

usage: python agreement.py DATA_TYPE SAMPLE CONTIG OUTPUT TOOL=PATH...
"""

import io
import subprocess
import sys
from itertools import combinations

import numpy as np
import pandas as pd
import pbzarr as pbz
import pyBigWig

data_type, sample, contig, output = sys.argv[1:5]
paths = dict(arg.split("=", 1) for arg in sys.argv[5:])


def load(path):
    if path.endswith(".d4"):
        # d4tools view, not pyd4: our pyd4 build misreads files that d4tools create writes.
        # view prints every interval of the contig, zeros included.
        view = subprocess.run(["d4tools", "view", path, contig], capture_output=True, check=True)
        bg = pd.read_csv(io.BytesIO(view.stdout), sep="\t", header=None, usecols=[1, 2, 3]).to_numpy()
        return np.repeat(bg[:, 2], bg[:, 1] - bg[:, 0]).astype(np.int64)
    if path.endswith(".bw"):
        with pyBigWig.open(path) as bw:
            values = bw.values(contig, 0, bw.chroms()[contig], numpy=True)
        return np.nan_to_num(values).astype(np.int64)  # no record = 0
    ds = pbz.open(path)["depth"].to_dataset(inherit=False)
    return ds.pbz.region(contig).values.astype(np.int64).ravel()


depth = {tool: load(path) for tool, path in paths.items()}

with open(output, "w") as out:
    out.write("data_type\tsample\tcontig\ttool_a\ttool_b\tcovered_bases\tdiffer_bases\tdiffer_pct\tmax_abs_diff\n")
    for a, b in combinations(depth, 2):
        diff = depth[a] - depth[b]
        covered = np.count_nonzero((depth[a] > 0) | (depth[b] > 0))
        differ = np.count_nonzero(diff)
        out.write(f"{data_type}\t{sample}\t{contig}\t{a}\t{b}\t{covered}\t{differ}"
                  f"\t{100 * differ / covered:.4f}\t{np.abs(diff).max()}\n")
