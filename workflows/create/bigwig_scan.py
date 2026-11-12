"""Exact pyBigWig means over full chromosomes or BED intervals.

Use the pyBigWig backend distributed with deepTools. Sum raw signal with
``stats(type="sum", exact=True)`` and divide by the requested interval length
so uncovered positions count as zero. No zoom summaries or per-base arrays.
"""

import argparse
import math
import os
from pathlib import Path
import sys

# Snakemake exports the whole rule's thread allocation to each xargs worker.
# Override it before pyBigWig imports NumPy and initializes a BLAS thread pool.
for variable in (
    "OMP_NUM_THREADS", "GOTO_NUM_THREADS", "OPENBLAS_NUM_THREADS",
    "OPENBLAS_DEFAULT_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[variable] = "1"

import pyBigWig


def bed_regions(filename, chrom_sizes):
    """Read BED in its original order, retaining duplicate intervals."""
    with Path(filename).open() as bed:
        for line_number, line in enumerate(bed, 1):
            fields = line.split()
            if not fields or fields[0].startswith("#") or fields[0] in {"track", "browser"}:
                continue
            location = f"{filename} line {line_number}"
            if len(fields) < 3 or not all(
                field.isascii() and field.isdecimal() for field in fields[1:3]
            ):
                raise ValueError(f"Invalid BED interval at {location}")
            chrom = fields[0]
            start, end = map(int, fields[1:3])
            if chrom not in chrom_sizes:
                raise ValueError(f"Unknown chromosome {chrom!r} at {location}")
            if not 0 <= start < end <= chrom_sizes[chrom]:
                raise ValueError(f"Invalid or out-of-range BED interval at {location}")
            yield chrom, start, end


def exact_mean(bw, chrom, start, end):
    """One exact sum query over the complete requested interval."""
    total = bw.stats(chrom, start, end, type="sum", exact=True)[0]
    if total is None:
        # pyBigWig also returns None for NaN signal and some corrupt data.
        # Check only ambiguous results: genuine uncovered regions have no
        # intervals, whereas corrupt reads raise and covered NaNs have records.
        if bw.intervals(chrom, start, end):
            raise ValueError(f"Missing or nonfinite BigWig sum for {chrom}:{start}-{end}")
        total = 0.0
    if not math.isfinite(total):
        raise ValueError(f"Nonfinite BigWig sum for {chrom}:{start}-{end}")
    return total / (end - start)


def scan(bigwig_path, bed_path=None, output=None):
    """Write headerless chromosome/start/end/mean rows to the output stream."""
    if output is None:
        output = sys.stdout
    with pyBigWig.open(str(bigwig_path)) as bw:
        if not bw.isBigWig():
            raise ValueError(f"Not a BigWig file: {bigwig_path}")
        chrom_sizes = bw.chroms()
        if bed_path is None:
            regions = ((chrom, 0, size) for chrom, size in chrom_sizes.items())
        else:
            regions = bed_regions(bed_path, chrom_sizes)
        for chrom, start, end in regions:
            if end <= start:
                raise ValueError(f"Zero-length chromosome {chrom!r}")
            mean = exact_mean(bw, chrom, start, end)
            print(f"{chrom}\t{start}\t{end}\t{mean:.17g}", file=output)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bigwig", metavar="BIGWIG")
    parser.add_argument("bed", metavar="BED", nargs="?")
    args = parser.parse_args(argv)
    try:
        scan(args.bigwig, args.bed)
        sys.stdout.flush()
    except (OSError, RuntimeError, ValueError) as error:
        print(f"bigwig_scan: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
