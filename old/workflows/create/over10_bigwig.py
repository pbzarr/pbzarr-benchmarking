"""Per-base fraction of samples with depth >= 10x (gnomAD's over_10), as a bigWig.

usage: python over10_bigwig.py OUT.bw WORKERS SAMPLE.bw [SAMPLE.bw ...]
"""

import sys
from functools import partial
from multiprocessing import Pool

import numpy as np
import pyBigWig


def over_10(paths, window):
    chrom, start, end = window
    covered = np.zeros(end - start)
    for path in paths:
        with pyBigWig.open(path) as bw:
            depth = np.nan_to_num(bw.values(chrom, start, end, numpy=True))  # no record = 0
        covered += depth >= 10
    return covered / len(paths)


if __name__ == "__main__":
    out, workers, paths = sys.argv[1], int(sys.argv[2]), sys.argv[3:]

    with pyBigWig.open(paths[0]) as bw:
        chroms = bw.chroms()
    windows = [(chrom, start, min(start + 1_000_000, length))
               for chrom, length in chroms.items()
               for start in range(0, length, 1_000_000)]

    bw = pyBigWig.open(out, "w")
    bw.addHeader(list(chroms.items()))
    with Pool(workers) as pool:
        # imap returns windows in order, which is the order bigWigs are written in
        for (chrom, start, end), values in zip(windows, pool.imap(partial(over_10, paths), windows)):
            bw.addEntries(chrom, start, values=values, span=1, step=1)
    bw.close()
