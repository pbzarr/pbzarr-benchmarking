#!/usr/bin/env python3
"""Time back-to-back trials of one command under a controlled page-cache state.

cold: before every trial, drop the page cache of each file under --data with
      posix_fadvise(DONTNEED) (Linux, no root; clean pages only).
warm: run the command once untimed, then time each trial back to back.
--before runs a shell command untimed ahead of every run (e.g. removing the
previous run's output).

Writes one TSV row per trial with Snakemake-benchmark-compatible fields:
s (wall seconds), max_rss (MiB, summed over the live process tree, sampled
every SAMPLE_SECONDS), io_in (MiB read from storage, from rusage block
inputs), cpu_time (user+sys seconds of all reaped descendants), and
mean_load (100 * cpu_time / s).
"""

import argparse
import csv
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

import psutil


SAMPLE_SECONDS = 0.05
FIELDS = ["trial", "cache", "s", "max_rss", "io_in", "cpu_time", "mean_load"]


def data_files(roots):
    for root in roots:
        if root.is_file():
            yield root
        for path in root.rglob("*"):
            if path.is_file() and not path.is_symlink():
                yield path


def evict(roots):
    os.sync()
    for path in data_files(roots):
        fd = os.open(path, os.O_RDONLY)
        try:
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
        finally:
            os.close(fd)


def tree_rss(process):
    total = 0
    try:
        members = [process, *process.children(recursive=True)]
    except psutil.NoSuchProcess:
        return 0
    for member in members:
        try:
            total += member.memory_info().rss
        except psutil.NoSuchProcess:
            pass
    return total


def run_once(command):
    """Run command to completion; return wall, peak tree RSS and rusage deltas."""
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.perf_counter()
    process = psutil.Popen(command)
    peak = 0
    while process.poll() is None:
        peak = max(peak, tree_rss(process))
        time.sleep(SAMPLE_SECONDS)
    wall = time.perf_counter() - start
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    if process.returncode != 0:
        raise subprocess.CalledProcessError(process.returncode, command)
    cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
    blocks = after.ru_inblock - before.ru_inblock
    return {
        "s": f"{wall:.4f}",
        "max_rss": f"{peak / 2**20:.2f}",
        "io_in": f"{blocks * 512 / 2**20:.2f}",
        "cpu_time": f"{cpu:.2f}",
        "mean_load": f"{100 * cpu / wall:.2f}" if wall > 0 else "0",
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, action="append", required=True,
                        help="File or directory whose cache state is controlled (repeatable)")
    parser.add_argument("--cache", choices=["cold", "warm"], required=True)
    parser.add_argument("--trials", type=int, default=5)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--before", help="Shell command run untimed before every run")
    parser.add_argument("command", nargs=argparse.REMAINDER,
                        help="Command to time, after --")
    args = parser.parse_args(argv)
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("missing command after --")
    if args.trials < 1:
        parser.error("--trials must be at least 1")
    if args.cache == "cold" and not hasattr(os, "posix_fadvise"):
        parser.error("cold cache needs os.posix_fadvise, which this platform lacks")

    def before():
        if args.before:
            subprocess.run(args.before, shell=True, check=True)

    try:
        rows = []
        if args.cache == "warm":
            before()
            run_once(command)
        for trial in range(1, args.trials + 1):
            before()
            if args.cache == "cold":
                evict(args.data)
            rows.append({"trial": trial, "cache": args.cache, **run_once(command)})
    except subprocess.CalledProcessError as error:
        print(f"timed_trials: command exited with status {error.returncode}", file=sys.stderr)
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
