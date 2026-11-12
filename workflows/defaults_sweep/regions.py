"""Write N random 10 kb regions on one contig as sorted BED (fixed seed).

usage: python regions.py CONTIG LENGTH N SEED > regions.bed
"""

import sys

import numpy as np

contig, length, n, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(seed)
starts = np.sort(rng.integers(0, length - 10_000, size=n))
for s in starts:
    print(f"{contig}\t{s}\t{s + 10_000}")
