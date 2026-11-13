#!/usr/bin/env python3
"""Collect trial benchmarks and du -sk results into two TSVs."""

import argparse
import csv
from pathlib import Path
import sys


TRIAL_FIELDS = ["dataset", "tool", "threads", "n_samples", "replicate"]
BENCHMARK_FIELDS = [
    "s", "h:m:s", "max_rss", "max_vms", "max_uss", "max_pss",
    "io_in", "io_out", "mean_load", "cpu_time",
]


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{path}: {len(rows)} rows")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "results", type=Path, nargs="?",
        default=Path(__file__).resolve().parents[2] / "results",
        help="Results root containing dataset directories (default: repo/results)",
    )
    results = parser.parse_args().results.resolve()
    if not results.is_dir():
        parser.error(f"Results directory does not exist: {results}")

    benchmarks, sizes = [], []
    metadata_fields = TRIAL_FIELDS + ["rule", "cache", "sample", "measurement", "source"]
    metric_fields = BENCHMARK_FIELDS.copy()

    for trial in sorted(results.glob("*/*/t*/n*/r*")):
        parts = trial.relative_to(results).parts
        if not trial.is_dir() or not all(part[1:].isdigit() for part in parts[2:]):
            continue
        metadata = dict(zip(TRIAL_FIELDS, parts))
        for key in ("threads", "n_samples", "replicate"):
            metadata[key] = int(metadata[key][1:])

        paths = [trial / "benchmark.tsv"]
        paths += sorted((trial / "benchmarks").rglob("*.tsv"))
        for path in paths:
            if not path.is_file():
                continue
            rule, sample, cache = "create", "", ""
            if path.parent != trial:
                relative = path.relative_to(trial / "benchmarks")
                rule = relative.stem if len(relative.parts) == 1 else relative.parent.as_posix()
                sample = "" if len(relative.parts) == 1 else relative.stem
            # Scan trials are named {analysis}.{cache}.tsv.
            rule, _, cache = rule.partition(".")
            with path.open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if not reader.fieldnames or "s" not in reader.fieldnames:
                    print(f"Skipping empty or invalid benchmark: {path}", file=sys.stderr)
                    continue
                for measurement, row in enumerate(reader, 1):
                    if None in row or any(value is None for value in row.values()):
                        print(f"Skipping incomplete row {measurement}: {path}", file=sys.stderr)
                        continue
                    # timed_trials.py rows carry their own trial number and
                    # cache mode; promote them to metadata.
                    measurement = int(row.pop("trial", measurement))
                    row_cache = row.pop("cache", cache)
                    if row_cache != cache:
                        print(f"Cache {row_cache!r} disagrees with filename: {path}", file=sys.stderr)
                        continue
                    # Preserve extended benchmark fields without replacing trial metadata.
                    metrics = {
                        f"benchmark_{key}" if key in metadata_fields else key: value
                        for key, value in row.items()
                    }
                    metric_fields.extend(key for key in metrics if key not in metric_fields)
                    benchmarks.append({
                        **metadata, "rule": rule, "cache": cache, "sample": sample,
                        "measurement": measurement, "source": str(path.relative_to(results)),
                        **metrics,
                    })

        disk = trial / "disk.tsv"
        if disk.is_file():
            for line in disk.read_text().splitlines():
                if not line.strip():
                    continue
                size, target = line.split(maxsplit=1)
                sizes.append({
                    **metadata, "size_kib": int(size), "size_gib": int(size) / 1024**2,
                    "target": target, "source": str(disk.relative_to(results)),
                })

    write_tsv(results / "benchmarks.tsv", benchmarks, metadata_fields + metric_fields)
    write_tsv(results / "sizes.tsv", sizes, TRIAL_FIELDS + ["size_kib", "size_gib", "target", "source"])


if __name__ == "__main__":
    main()
