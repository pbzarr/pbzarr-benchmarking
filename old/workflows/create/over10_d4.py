"""Per-base number of samples with depth >= 10x (gnomAD's over_10 times n), as a D4 file.

usage: python over10_d4.py COHORT.d4 OUT.d4 WORKERS
"""

import sys
from functools import partial
from multiprocessing import Pool

import pyd4


def over_10(path, window):
    depth = pyd4.D4File(path).open_all_tracks()[window]  # samples x positions
    return (depth >= 10).sum(axis=0)  # D4 stores integers, so keep the count


if __name__ == "__main__":
    path, out, workers = sys.argv[1], sys.argv[2], int(sys.argv[3])

    chroms = pyd4.D4File(path).chroms()
    windows = [(chrom, start, min(start + 1_000_000, length))
               for chrom, length in chroms
               for start in range(0, length, 1_000_000)]

    writer = pyd4.D4Builder(out).add_chroms(chroms).set_dict_bits(8).get_writer()
    with Pool(workers) as pool:
        for (chrom, start, end), counts in zip(windows, pool.imap(partial(over_10, path), windows)):
            writer.write_np_array(chrom, start, counts)
    writer.close()
