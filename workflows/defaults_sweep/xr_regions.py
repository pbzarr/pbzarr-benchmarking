"""Read each BED region of a one-contig PBZ store with plain xarray (zarr-python).

mode all: every sample; mode one: the first sample only.
usage: python xr_regions.py STORE REGIONS.bed all|one
"""

import sys

import xarray as xr

store, bed, mode = sys.argv[1], sys.argv[2], sys.argv[3]
ds = xr.open_zarr(store, chunks=None, consolidated=True)
depth = ds["depth"]
total = 0
for line in open(bed):
    _, s, e = line.split()[:3]
    sel = depth[int(s):int(e)] if mode == "all" else depth[int(s):int(e), 0]
    total += int(sel.values.sum())
print(total)
