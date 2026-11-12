"""Print one TSV row of storage facts for STORE/depth.

Columns: files, bytes, bytes_per_cell, inner_chunks, inner_stored,
empty_frac, shards, mean_shard_bytes, index_bytes.
Empty inner chunks are found from the shard index (offset and size both
2^64-1) or, without sharding, from missing chunk files.
usage: python measure_store.py STORE [--header]
"""

import json
import math
import os
import sys

import numpy as np

EMPTY = np.uint64(2**64 - 1)


def main():
    store = sys.argv[1]
    if "--header" in sys.argv:
        print("files\tbytes\tbytes_per_cell\tinner_chunks\tinner_stored\tempty_frac\tshards\tmean_shard_bytes\tindex_bytes")
        return
    meta = json.load(open(os.path.join(store, "depth", "zarr.json")))
    shape = meta["shape"]
    outer = meta["chunk_grid"]["configuration"]["chunk_shape"]
    shard = next((c for c in meta["codecs"] if c["name"] == "sharding_indexed"), None)
    files = []
    for root, _, names in os.walk(os.path.join(store, "depth", "c")):
        files += [os.path.join(root, f) for f in names]
    sizes = [os.path.getsize(f) for f in files]
    total = sum(sizes)
    cells = shape[0] * shape[1]
    if shard is None:
        inner = outer
        inner_chunks = math.ceil(shape[0] / inner[0]) * math.ceil(shape[1] / inner[1])
        stored = len(files)
        index_bytes = 0
    else:
        inner = shard["configuration"]["chunk_shape"]
        per_shard = [outer[d] // inner[d] for d in range(2)]
        entries = per_shard[0] * per_shard[1]
        has_crc = any(c["name"] == "crc32c" for c in shard["configuration"]["index_codecs"])
        index_bytes = entries * 16 + (4 if has_crc else 0)
        inner_chunks = math.ceil(shape[0] / inner[0]) * math.ceil(shape[1] / inner[1])
        stored = 0
        for f in files:
            with open(f, "rb") as fh:
                fh.seek(-index_bytes, os.SEEK_END)
                idx = np.frombuffer(fh.read(entries * 16), dtype="<u8").reshape(-1, 2)
            stored += int(np.count_nonzero(~((idx[:, 0] == EMPTY) & (idx[:, 1] == EMPTY))))
        # entries past the array end are always empty; count only real chunks
        stored = min(stored, inner_chunks)
    empty = 1 - stored / inner_chunks
    mean_shard = total / len(files) if files else 0
    print(f"{len(files)}\t{total}\t{total / cells:.5f}\t{inner_chunks}\t{stored}\t{empty:.4f}\t"
          f"{len(files) if shard else 0}\t{mean_shard:.0f}\t{index_bytes}")


if __name__ == "__main__":
    main()
