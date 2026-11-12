"""over10_pbz.py without dask: a plain loop over shard-aligned blocks.

usage: python over10_pbz_loop.py COHORT.pbz OUT.pbz THREADS
"""

import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pbzarr as pbz

if __name__ == "__main__":
    path, out, threads = sys.argv[1], sys.argv[2], int(sys.argv[3])

    ds = pbz.open(path)
    depth = ds.values("depth")  # positions x samples
    step = depth.shards[0]
    over_10 = pbz.create(out, reference=ds.reference).empty("over_10", dtype="float32", chunks=step)

    def block(start):
        # numpy and the zarr decoders release the GIL, so threads run in parallel
        over_10[start:start + step] = (depth[start:start + step] >= 10).mean(axis=1, dtype=np.float32)

    with ThreadPoolExecutor(threads) as pool:
        list(pool.map(block, range(0, depth.shape[0], step)))
