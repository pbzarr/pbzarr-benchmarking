"""Per-base coverage summary (gnomAD style) as D4 files: total depth (gnomAD's
total_DP) and the number of samples with depth >= each threshold, one D4 file per
output. D4 stores integers, so it keeps sums and counts rather than the mean and
fractions.

Sized to N like workflows/cohort/over10_d4.py: the window is min(1 Mb, 128 Mb / N)
bases, the count files have enough dictionary bits for 0..N, and total_DP has 16
bits (totals up to 65,535, e.g. N = 1,024 at 64x, stay in the primary table).

usage: python summary_d4.py COHORT.d4 OUT_DIR WORKERS
"""

import sys
from functools import partial
from math import ceil, log2
from multiprocessing import Pool
from pathlib import Path

import pyd4

THRESHOLDS = [1, 5, 10, 15, 20, 25, 30, 50, 100]
NAMES = ["total_dp"] + [f"over_{t}" for t in THRESHOLDS]


def summarize(path, window):
    depth = pyd4.D4File(path).open_all_tracks()[window]  # samples x positions
    return [depth.sum(axis=0), *((depth >= t).sum(axis=0) for t in THRESHOLDS)]


if __name__ == "__main__":
    path, out, workers = sys.argv[1], Path(sys.argv[2]), int(sys.argv[3])

    chroms = pyd4.D4File(path).chroms()
    n = len(pyd4.D4File(path).list_tracks())  # number of samples in the file
    size = min(1_000_000, 128_000_000 // n)
    windows = [(chrom, start, min(start + size, length))
               for chrom, length in chroms
               for start in range(0, length, size)]

    out.mkdir()
    bits = [16] + [ceil(log2(n + 1))] * len(THRESHOLDS)
    writers = [pyd4.D4Builder(str(out / f"{name}.d4")).add_chroms(chroms).set_dict_bits(b).get_writer()
               for name, b in zip(NAMES, bits)]
    with Pool(workers) as pool:
        for (chrom, start, end), tracks in zip(windows, pool.imap(partial(summarize, path), windows)):
            for writer, values in zip(writers, tracks):
                writer.write_np_array(chrom, start, values)
    for writer in writers:
        writer.close()
