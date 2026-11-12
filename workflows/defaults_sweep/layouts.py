"""Print the layouts of the one-hour defaults sweep as TSV.

Columns: name, grid, chunk, shard, a2a, a2b, b2b. chunk and shard are
comma-separated shapes for zarrs_reencode -c and -s (shard "-" = no
sharding). a2a, a2b and b2b are codec JSON ("-" = default).
usage: python layouts.py N_SAMPLES > layouts.tsv
"""

import json
import sys

S = 16  # samples per shard
BYTES = {"name": "bytes", "configuration": {"endian": "little"}}


def blosc(cname, level, blocksize=0):
    return [{"name": "blosc", "configuration": {
        "cname": cname, "clevel": level, "shuffle": "shuffle", "typesize": 4, "blocksize": blocksize}}]


CODECS = {
    "bzstd1": ("", BYTES, blosc("zstd", 1)),
    "bzstd3": ("", BYTES, blosc("zstd", 3)),
    "bzstd5": ("", BYTES, blosc("zstd", 5)),
    "zstd3": ("", BYTES, [{"name": "zstd", "configuration": {"level": 3, "checksum": False}}]),
    "blz4": ("", BYTES, blosc("lz4", 5)),
    "pcodec8": ("", {"name": "numcodecs.pcodec", "configuration": {"level": 8}}, []),
}


def row(name, grid, chunk, shard, a2a, a2b, b2b):
    enc = lambda v: json.dumps(v, separators=(",", ":")) if v not in ("", []) else "-"
    return "\t".join([name, grid, chunk, shard or "-", enc(a2a), enc(a2b), enc(b2b)])


def sharded(pc, cps):
    return f"{pc},1", f"{pc * cps},{S}"


def main():
    n = int(sys.argv[1])
    rows = []
    # codec grid at 65,536 positions per chunk and 256 chunks per shard
    for key, (a2a, a2b, b2b) in CODECS.items():
        rows.append(row(f"codec-{key}", "codec", *sharded(65536, 256), a2a, a2b, b2b))
    rows.append(row("codec-bzstd3-wholeblock", "codec", *sharded(65536, 256), "", BYTES,
                    blosc("zstd", 3, blocksize=65536 * 4)))
    rows.append(row("ref-control", "codec", f"1000000,{min(16, n)}", "", "", BYTES, blosc("zstd", 5)))
    rows.append(row("ref-bench", "codec", "16384,8", "1048576,8",
                    [{"name": "transpose", "configuration": {"order": [1, 0]}}],
                    {"name": "numcodecs.pcodec", "configuration": {"level": 8}}, []))
    # geometry grid for two codecs
    for key in ("bzstd3", "pcodec8"):
        a2a, a2b, b2b = CODECS[key]
        for pc in (16384, 262144, 1000000):
            rows.append(row(f"geo-{key}-pc{pc}-cps256", "geometry", *sharded(pc, 256), a2a, a2b, b2b))
        # at 16,384 positions a 4,096-chunk shard still fits in chr20, so the index size is real
        for cps in (64, 1024, 4096):
            rows.append(row(f"geo-{key}-pc16384-cps{cps}", "geometry", *sharded(16384, cps), a2a, a2b, b2b))
    print("\n".join(rows))


if __name__ == "__main__":
    main()
