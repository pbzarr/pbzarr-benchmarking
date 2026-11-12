"""Per-base coverage summary (gnomAD style) as bigWigs: mean depth and the
fraction of samples with depth >= each threshold, one bigWig per output.

usage: python summary_bigwig.py OUT_DIR WORKERS SAMPLE.bw [SAMPLE.bw ...]
"""

import sys
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pyBigWig

THRESHOLDS = [1, 5, 10, 15, 20, 25, 30, 50, 100]
NAMES = ["mean"] + [f"over_{t}" for t in THRESHOLDS]


def summarize(paths, window):
    chrom, start, end = window
    total = np.zeros(end - start)
    over = np.zeros((len(THRESHOLDS), end - start))
    for path in paths:
        with pyBigWig.open(path) as bw:
            depth = np.nan_to_num(bw.values(chrom, start, end, numpy=True))  # no record = 0
        total += depth
        for i, t in enumerate(THRESHOLDS):
            over[i] += depth >= t
    return [total / len(paths), *(over / len(paths))]


if __name__ == "__main__":
    out, workers, paths = Path(sys.argv[1]), int(sys.argv[2]), sys.argv[3:]

    with pyBigWig.open(paths[0]) as bw:
        chroms = bw.chroms()
    windows = [(chrom, start, min(start + 1_000_000, length))
               for chrom, length in chroms.items()
               for start in range(0, length, 1_000_000)]

    out.mkdir()
    writers = [pyBigWig.open(str(out / f"{name}.bw"), "w") for name in NAMES]
    for bw in writers:
        bw.addHeader(list(chroms.items()))
    with Pool(workers) as pool:
        # imap returns windows in order, which is the order bigWigs are written in
        for (chrom, start, end), tracks in zip(windows, pool.imap(partial(summarize, paths), windows)):
            for bw, values in zip(writers, tracks):
                bw.addEntries(chrom, start, values=values, span=1, step=1)
    for bw in writers:
        bw.close()
