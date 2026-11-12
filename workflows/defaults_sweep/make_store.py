"""Make a 0.5 PBZ store around a depth array written by zarrs_reencode.

pre:  python make_store.py pre META_JSON STORE
      creates the store (reference coordinate) and the sample label array.
post: python make_store.py post STORE
      writes consolidated metadata after STORE/depth exists.
Run with the Python env that has pbzarr 0.8.
"""

import json
import sys

import numpy as np
import pbzarr
import zarr


def pre(meta_path, store):
    meta = json.loads(open(meta_path).read())
    ds = pbzarr.create(store, reference={meta["contig"]: meta["length"]})
    ds.close()
    labels = np.array(meta["samples"], dtype=str)
    arr = zarr.create_array(store, name="sample", shape=labels.shape, chunks=labels.shape,
                            dtype=str, dimension_names=["sample"])
    arr[:] = labels


def post(store):
    zarr.consolidate_metadata(store)


if __name__ == "__main__":
    if sys.argv[1] == "pre":
        pre(sys.argv[2], sys.argv[3])
    else:
        post(sys.argv[2])
