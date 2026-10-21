"""How many D4 or BigWig jobs to run at once on the node, and with how many threads.

  packing.py choose RUNS.tsv [--cores 256] [--mem-gb 1400]
      Read part 1's runs.tsv and print the packing as YAML, keyed d4 (the
      mosdepth rows) and bigwig (the Megadepth rows).
  packing.py for PACKING.yaml TOOL N
      Print "THREADS JOBS" for a cohort of N samples.

Per thread count, part 1 gives medians over samples and repeats of wall_s,
cpu_s (user_s + sys_s) and max RSS. From those:
  fastest_threads   the thread count with the lowest wall_s.
  efficient_threads the lowest cpu_s among counts with wall_s at most twice the
                    fastest.
  efficient_jobs    min(cores / busy, memory budget / max RSS), rounded down,
                    where busy = cpu_s / wall_s at the efficient thread count.
For a cohort of N samples: if N >= efficient_jobs, run efficient_jobs jobs of
efficient_threads threads. Otherwise run N jobs of fastest_threads threads.
"""

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict

import yaml

# Tool name in the packing file -> tool name in part 1's runs.tsv.
RUN_TOOLS = {"d4": "mosdepth", "bigwig": "megadepth"}


def medians_by_threads(rows):
    """{threads: (wall_s, cpu_s, max_rss_kb)}, medians over samples and repeats."""
    groups = defaultdict(list)
    for row in rows:
        groups[int(row["threads"])].append(row)
    return {
        threads: (
            statistics.median(float(row["wall_s"]) for row in group),
            statistics.median(float(row["user_s"]) + float(row["sys_s"]) for row in group),
            statistics.median(float(row["max_rss_kb"]) for row in group),
        )
        for threads, group in groups.items()
    }


def choose_packing(medians, cores, mem_gb):
    wall = {threads: m[0] for threads, m in medians.items()}
    cpu = {threads: m[1] for threads, m in medians.items()}
    fastest = min(wall, key=lambda threads: (wall[threads], threads))
    close = [threads for threads in wall if wall[threads] <= 2 * wall[fastest]]
    efficient = min(close, key=lambda threads: (cpu[threads], threads))
    busy = cpu[efficient] / wall[efficient]
    mem_kb = mem_gb * 1024 * 1024
    jobs = min(math.floor(cores / busy), math.floor(mem_kb / medians[efficient][2]))
    return {"efficient_threads": efficient, "efficient_jobs": jobs, "fastest_threads": fastest}


def choose(args):
    with open(args.runs) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    packing = {}
    for tool, run_tool in RUN_TOOLS.items():
        tool_rows = [row for row in rows if row["tool"] == run_tool]
        if not tool_rows:
            sys.exit(f"packing.py: no {run_tool} rows in {args.runs}")
        packing[tool] = choose_packing(medians_by_threads(tool_rows), args.cores, args.mem_gb)
    print(yaml.safe_dump(packing, sort_keys=False), end="")


def for_cohort(args):
    with open(args.packing) as handle:
        entry = yaml.safe_load(handle)[args.tool]
    if args.n >= entry["efficient_jobs"]:
        print(entry["efficient_threads"], entry["efficient_jobs"])
    else:
        print(entry["fastest_threads"], args.n)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)

    choose_parser = commands.add_parser("choose", help="pick the packing from part 1's runs.tsv")
    choose_parser.add_argument("runs", help="part 1's runs.tsv")
    choose_parser.add_argument("--cores", type=int, default=256)
    choose_parser.add_argument("--mem-gb", type=float, default=1400,
                               help="memory for the jobs, in GiB")
    choose_parser.set_defaults(run=choose)

    for_parser = commands.add_parser("for", help="print THREADS JOBS for a cohort of N samples")
    for_parser.add_argument("packing", help="YAML written by `choose`")
    for_parser.add_argument("tool", choices=sorted(RUN_TOOLS))
    for_parser.add_argument("n", type=int)
    for_parser.set_defaults(run=for_cohort)

    args = parser.parse_args()
    args.run(args)


if __name__ == "__main__":
    main()
