"""Copy one contig of a 0.4 PBZ depth group to a plain uncompressed zarr array.

The copy has chunks [65536, 1], so zarrs_reencode reads little extra for any
target layout. Also writes meta.json (contig, length, samples).
usage: python extract_contig.py SRC_DEPTH_GROUP CONTIG OUT_DIR WORKERS [ROWS]
ROWS limits the copy to the first ROWS positions (for a smoke test).
"""

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import zarr

BLOCK = 1 << 20
CHUNK = 65536


def copy_block(args):
    src, out, start, a, b = args
    values = zarr.open_array(f"{src}/values", mode="r")
    raw = zarr.open_array(out, mode="r+")
    raw[a:b, :] = values[start + a:start + b, :]
    return b - a


def main():
    src, contig, out_dir, workers = sys.argv[1], sys.argv[2], Path(sys.argv[3]), int(sys.argv[4])
    rows = int(sys.argv[5]) if len(sys.argv) > 5 else None
    g = zarr.open_group(src, mode="r")
    names = [str(x) for x in g["contigs"][:]]
    offsets = g["offsets"][:]
    samples = [str(x) for x in g["sample"][:]]
    i = names.index(contig)
    start, length = int(offsets[i]), int(offsets[i + 1] - offsets[i])
    if rows:
        length = min(length, rows)
    n = len(samples)
    out_dir.mkdir(parents=True, exist_ok=True)
    out = str(out_dir / "raw")
    zarr.create_array(out, shape=(length, n), chunks=(CHUNK, 1), dtype="int32",
                      fill_value=0, compressors=None, filters=None, overwrite=True)
    blocks = [(src, out, start, a, min(a + BLOCK, length)) for a in range(0, length, BLOCK)]
    with Pool(workers) as pool:
        done = sum(pool.imap_unordered(copy_block, blocks))
    assert done == length
    (out_dir / "meta.json").write_text(json.dumps(
        {"contig": contig, "length": length, "samples": samples, "source": src, "start": start}))
    print(f"{contig}: {length} x {n} copied")


if __name__ == "__main__":
    main()
