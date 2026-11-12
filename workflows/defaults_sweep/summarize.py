"""Summarize a defaults sweep run as Markdown tables, with the decision rules.

usage: python summarize.py RESULTS_DIR > summary.md
"""

import csv
import math
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

GENOME = 3_217_346_917  # GRCh38 analysis set with decoys, as in the dense input
S = 16


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if path.exists() else []


def med(xs):
    return statistics.median(xs) if xs else float("nan")


def _unused_geometry(name):
    m = re.search(r"pc(\d+)-cps(\d+)", name)
    if m:
        return int(m.group(1)), int(m.group(2))
    if name.startswith("codec-"):
        return 65536, 256
    return None


def main():
    out = Path(sys.argv[1])
    enc = {(r["input"], r["name"]): r for r in tsv(out / "encode.tsv")}

    dec = defaultdict(list)
    for line in open(out / "decode.tsv") if (out / "decode.tsv").exists() else []:
        st, _, text = line.rstrip("\n").split("\t", 2)
        m = re.search(r"@ ([\d.]+)GB/s", text)
        if m:
            p = Path(st)
            dec[(p.parent.name, p.stem)].append(float(m.group(1)) / 4)

    reads = defaultdict(list)
    for r in tsv(out / "reads.tsv"):
        reads[(r["location"], r["input"], r["name"], r["test"])].append(float(r["s"]))

    small = {}
    ctx = (out / "context.txt").read_text() if (out / "context.txt").exists() else ""
    m = re.search(r"small_cells: (.*)", ctx)
    if m:
        parts = m.group(1).split()
        small = {parts[k]: int(parts[k + 1]) for k in range(0, len(parts) - 1, 2)}
    spd = defaultdict(list)
    for r in tsv(out / "encode_speed.tsv"):
        spd[(r["input"], r["name"])].append(float(r["wall_s"]))

    def enc_speed(r):
        """MB of raw int32 encoded per second per thread (4 pinned threads)."""
        walls = spd[(r["input"], r["name"])]
        if not walls or r["input"] not in small:
            return float("nan")
        return small[r["input"]] * 4 / 1e6 / med(walls) / 4

    print("# Defaults sweep summary\n")
    for inp in sorted({k[0] for k in enc}):
        rows = [r for (i, _), r in sorted(enc.items()) if i == inp]
        # ref-today is today's [1M x 16] unsharded layout (built by control.sh);
        # ref-control kept the input chunks [65536 x 1] because zarrs_reencode ignores -c without -s
        ctrl = enc.get((inp, "ref-today")) or enc.get((inp, "ref-control"))
        print(f"## {inp}: storage, encode and decode\n")
        print("| layout | status | B/cell | vs control | files | empty | index B | mean file MB | encode MB/s/thread | decode GB/s/thread |")
        print("|---|---|---|---|---|---|---|---|---|---|")
        for r in rows:
            if r["status"] != "ok":
                print(f"| {r['name']} | {r['status']} | | | | | | | | |")
                continue
            ratio = float(ctrl["bytes"]) / float(r["bytes"]) if ctrl and ctrl["status"] == "ok" else float("nan")
            print(f"| {r['name']} | ok | {float(r['bytes_per_cell']):.4f} | {ratio:.2f}x | {r['files']} | "
                  f"{float(r['empty_frac']):.3f} | {r['index_bytes']} | {float(r['mean_shard_bytes']) / 1e6:.2f} | "
                  f"{enc_speed(r):.0f} | {med(dec[(inp, r['name'])]):.2f} |")
        print()

        # codec rule, against today's layout and against today's codec on the new layout
        for base in ("ref-today", "codec-bzstd5"):
            b = enc.get((inp, base))
            if not b or b["status"] != "ok":
                continue
            c_dec, c_enc = med(dec[(inp, base)]), enc_speed(b)
            enc_ok = (lambda r: True) if math.isnan(c_enc) else (lambda r: enc_speed(r) >= c_enc / 2)
            ok = [r for r in rows if r["status"] == "ok" and r["name"].startswith("codec-")
                  and med(dec[(inp, r["name"])]) >= c_dec and enc_ok(r)]
            pick = min(ok, key=lambda r: float(r["bytes_per_cell"]), default=None)
            print(f"Codec rule against {base} (decode >= {c_dec:.2f} GB/s/thread, encode >= {c_enc / 2:.0f} MB/s/thread): "
                  f"**{pick['name'] if pick else 'none passes'}**\n")

    print("## Reads (median seconds)\n")
    tests = ["contig_one", "reg_all", "reg_one", "xr_all", "xr_one"]
    print("| location | input | layout | " + " | ".join(tests) + " |")
    print("|---|---|---|" + "---|" * len(tests))
    for loc, inp, name in sorted({k[:3] for k in reads}):
        cells = [f"{med(reads[(loc, inp, name, t)]):.2f}" if reads[(loc, inp, name, t)] else "" for t in tests]
        print(f"| {loc} | {inp} | {name} | " + " | ".join(cells) + " |")
    print()

    print("## Geometry rules (Lustre reads)\n")
    for inp in sorted({k[1] for k in reads}):
        for codec in ("bzstd3", "pcodec8"):
            def t(name, test):
                return med(reads[("lustre", inp, name, test)])
            base_name = "codec-" + codec
            pcs = {16384: f"geo-{codec}-pc16384-cps256", 65536: base_name,
                   262144: f"geo-{codec}-pc262144-cps256", 1000000: f"geo-{codec}-pc1000000-cps256"}
            ref = t(pcs[1000000], "contig_one")
            good = [pc for pc, n in sorted(pcs.items()) if t(n, "contig_one") <= ref * 1.10]
            cpss = {64: f"geo-{codec}-pc16384-cps64", 256: f"geo-{codec}-pc16384-cps256",
                    1024: f"geo-{codec}-pc16384-cps1024", 4096: f"geo-{codec}-pc16384-cps4096"}
            base = t(cpss[256], "reg_all")
            okc = [c for c, n in sorted(cpss.items()) if t(n, "reg_all") <= base * 1.10]
            print(f"- {inp}, {codec}: smallest pos_chunk with whole-contig stat within 10% of 1M: "
                  f"**{good[0] if good else 'none'}**; largest chunks per shard with region read within 10% of 256: "
                  f"**{okc[-1] if okc else 'none'}**")
    print()

    ft = tsv(out / "filetree.tsv")
    if ft:
        print("## File operations on Lustre (median seconds; ms per file)\n")
        ops = ["create", "find", "du", "cp", "rm"]
        per_file = {}
        print("| files | " + " | ".join(ops) + " |")
        print("|---|" + "---|" * len(ops))
        for n in sorted({int(r["files"]) for r in ft}):
            vals = []
            for op in ops:
                xs = [float(r["seconds"]) for r in ft if int(r["files"]) == n and r["op"] == op]
                vals.append(f"{med(xs):.2f} ({1000 * med(xs) / n:.3f})" if xs else "")
                if xs:
                    per_file[op] = 1000 * med(xs) / n
            print(f"| {n} | " + " | ".join(vals) + " |")
        print("\nWhole-genome files and estimated time (largest tree's ms per file):\n")
        print("| pos_chunk | chunks/shard | N | files | find s | du s | cp s | rm s |")
        print("|---|---|---|---|---|---|---|---|")
        for pc, cps in [(16384, 64), (16384, 256), (16384, 1024), (16384, 4096), (65536, 256), (65536, 1024)]:
            for n in (128, 1024, 3200):
                files = math.ceil(GENOME / (pc * cps)) * math.ceil(n / S)
                est = [f"{files * per_file.get(op, float('nan')) / 1000:.1f}" for op in ("find", "du", "cp", "rm")]
                print(f"| {pc} | {cps} | {n} | {files:,} | " + " | ".join(est) + " |")
        print("\nControl layout today ([1M x 16], no shard): files = 3,218 x ceil(N/16).")


if __name__ == "__main__":
    main()
