"""Time file operations on a tree of N small files laid out like a Zarr store.

Files go in DIR/c/<row>/<col>, 64 per row directory (the column shards of
1,024 samples at 16 samples per shard), 1 KiB each. Prints TSV rows:
files, trial, op, seconds. Ops: create, find, du, cp, rm_copy, rm (no cp
above 200,000 files, to keep the run short).
usage: python filetree.py DIR N TRIALS [--header]
"""

import os
import shutil
import subprocess
import sys
import time

COLS = 64
PAYLOAD = os.urandom(1024)


def timed(op, fn, rows, n, trial):
    t = time.perf_counter()
    fn()
    rows.append(f"{n}\t{trial}\t{op}\t{time.perf_counter() - t:.3f}")


def create(root, n):
    for i in range(n):
        d = os.path.join(root, "c", str(i // COLS))
        if i % COLS == 0:
            os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, str(i % COLS)), "wb") as fh:
            fh.write(PAYLOAD)


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL)


def main():
    if "--header" in sys.argv:
        print("files\ttrial\top\tseconds")
        return
    base, n, trials = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    rows = []
    for trial in range(1, trials + 1):
        root = os.path.join(base, f"t{n}-{trial}")
        copy = root + "-copy"
        for p in (root, copy):
            shutil.rmtree(p, ignore_errors=True)
        timed("create", lambda: create(root, n), rows, n, trial)
        timed("find", lambda: sh(f"find {root} -type f | wc -l"), rows, n, trial)
        timed("du", lambda: sh(f"du -s {root}"), rows, n, trial)
        if n <= 200_000:
            timed("cp", lambda: sh(f"cp -r {root} {copy}"), rows, n, trial)
            timed("rm_copy", lambda: sh(f"rm -rf {copy}"), rows, n, trial)
        timed("rm", lambda: sh(f"rm -rf {root}"), rows, n, trial)
        print("\n".join(rows), flush=True)
        rows.clear()


if __name__ == "__main__":
    main()
