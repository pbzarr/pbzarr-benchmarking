"""Per-base number of samples with depth >= 10x (gnomAD's over_10 times n), as a D4 file.

Copy of workflows/create/over10_d4.py for the cohort benchmark, with two changes for
N up to 1024: the window shrinks to min(1 Mb, 128 Mb / N) bases, so a worker holds
about the memory it holds at N = 128, and the dictionary has enough bits for counts 0..N.

usage: python over10_d4.py COHORT.d4 OUT.d4 WORKERS
"""

import sys
from functools import partial
from math import ceil, log2
from multiprocessing import Pool

import pyd4


def over_10(path, window):
    depth = pyd4.D4File(path).open_all_tracks()[window]  # samples x positions
    return (depth >= 10).sum(axis=0)  # D4 stores integers, so keep the count


if __name__ == "__main__":
    path, out, workers = sys.argv[1], sys.argv[2], int(sys.argv[3])

    chroms = pyd4.D4File(path).chroms()
    n = len(pyd4.D4File(path).list_tracks())  # number of samples in the file
    size = min(1_000_000, 128_000_000 // n)
    windows = [(chrom, start, min(start + size, length))
               for chrom, length in chroms
               for start in range(0, length, size)]

    writer = pyd4.D4Builder(out).add_chroms(chroms).set_dict_bits(ceil(log2(n + 1))).get_writer()
    with Pool(workers) as pool:
        for (chrom, start, end), counts in zip(windows, pool.imap(partial(over_10, path), windows)):
            writer.write_np_array(chrom, start, counts)
    writer.close()
