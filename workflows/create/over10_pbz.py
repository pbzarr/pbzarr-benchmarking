"""Per-base fraction of samples with depth >= 10x (gnomAD's over_10), as a PBZ track.

usage: python over10_pbz.py COHORT.pbz OUT.pbz WORKERS
"""

import sys

import pbzarr as pbz
from pbzarr.cluster import DaskClient

if __name__ == "__main__":
    path, out, workers = sys.argv[1], sys.argv[2], int(sys.argv[3])

    # processes sidestep the GIL in zarr's reads; a few threads each share the sample sums
    with DaskClient(workers=max(1, workers // 4), threads=4).get_client():
        ds = pbz.open(path)
        over_10 = (ds["depth"] >= 10).mean("sample").astype("float32")  # lazy, per base

        out = pbz.create(out, reference=ds.reference)
        out.write({"over_10": over_10})
