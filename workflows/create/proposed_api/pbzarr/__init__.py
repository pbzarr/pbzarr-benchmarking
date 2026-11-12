"""The proposed pbzarr API (open, ds[track], create, write) on format 0.4 stores.

A stand-in until the proposed Dataset/Track API ships. It runs the same lazy
xarray/dask operations: tracks open as dask arrays chunked by shard, and
writes copy the source's genome coordinates and PBZ metadata so the result is
a valid format 0.4 track. Put this directory on PYTHONPATH to use it as
`import pbzarr`.
"""

import dask
import numpy as np
import xarray as xr
import zarr


class Dataset:
    def __init__(self, path, reference):
        self.path = path
        self.reference = reference

    def __getitem__(self, track):
        values = xr.open_zarr(self.path, group=track, consolidated=False, chunks=None)["values"]
        shards = values.encoding.get("shards") or values.encoding["chunks"]
        chunks = dict(zip(values.dims, shards))
        # Half a shard of positions per dask chunk: a full int32 shard is 32 MiB, glibc's
        # largest mmap threshold, so every chunk would be a fresh mmap and page faults
        # dominate under many workers. Run with GLIBC_TUNABLES (see the over10 rules).
        chunks["position"] = max(1, chunks["position"] // 2)
        return values.chunk(chunks)

    def values(self, track):
        """A track's values as a plain zarr array: numpy in, numpy out, no dask."""
        return zarr.open_array(self.path, path=f"{track}/values")

    def empty(self, name, *, dtype, chunks):
        """Write a track's genome coordinates and metadata; return its values to fill."""
        genome = self.reference["genome"].compute()
        genome.to_zarr(self.path, group=name, consolidated=False)
        group = zarr.open_group(self.path, path=name)
        group.attrs.update(self.reference["track_attrs"])
        return group.create_array("values", shape=(int(genome["offsets"][-1]),), dtype=dtype,
                                  chunks=(chunks,), fill_value=np.nan,
                                  dimension_names=["position"])

    def write(self, tracks):
        # One compute for all tracks, so each input chunk is read once.
        writes = []
        for name, values in tracks.items():
            track = self.reference["genome"].assign(values=values.reset_coords(drop=True))
            track.attrs = self.reference["track_attrs"]
            writes.append(track.to_zarr(self.path, group=name, consolidated=False, compute=False))
        dask.compute(*writes)


def open(path):
    root = zarr.open_group(path, mode="r")
    first = next(name for name, _ in root.groups())
    track = xr.open_zarr(path, group=first, consolidated=False)
    reference = {
        "genome": track[["contigs", "offsets"]].reset_coords(),
        "track_attrs": dict(track.attrs),
        "store_attrs": dict(root.attrs),
    }
    return Dataset(path, reference)


def create(path, *, reference):
    zarr.open_group(path, mode="w-", attributes=reference["store_attrs"])
    return Dataset(path, reference)
