"""Per-base coverage summary (gnomAD style) as PBZ tracks: mean depth and the
fraction of samples with depth >= each threshold.

usage: python summary_pbz.py COHORT.pbz OUT.pbz WORKERS
"""

import sys

import pbzarr as pbz
from pbzarr.cluster import DaskClient

THRESHOLDS = [1, 5, 10, 15, 20, 25, 30, 50, 100]

if __name__ == "__main__":
    path, out, workers = sys.argv[1], sys.argv[2], int(sys.argv[3])

    # processes sidestep the GIL in zarr's reads; a few threads each share the sample sums
    with DaskClient(workers=max(1, workers // 4), threads=4).get_client():
        ds = pbz.open(path)
        depth = ds["depth"]
        summary = {"mean": depth.mean("sample").astype("float32")}
        for t in THRESHOLDS:
            summary[f"over_{t}"] = (depth >= t).mean("sample").astype("float32")

        out = pbz.create(out, reference=ds.reference)
        out.write(summary)  # all tracks in one pass over the input
