# Shared by the three cohort workflows (pbz.smk, d4.smk, bigwig.smk). Only the
# cohort workflows read this file.
#
# The batch driver (part2.sbatch) runs one workflow per tool in a working
# directory on the node's local disk (snakemake --directory). It passes these
# values with --config, on top of the config file:
#   n_samples, samples_file, results (this job's results folder), build,
#   reference_local (the node-local copy of the reference), and
#   threads_per_job (D4 and BigWig only).
#
# What the driver must also pass to snakemake, because Snakemake schedules by
# the threads and resources the rules declare:
#   --cores      at least max(CORES, JOBS * threads_per_job). Each D4 or BigWig
#                job declares threads_per_job, so with fewer cores fewer than
#                JOBS jobs run at once. (The guard below checks the CORES part.)
#   --resources  pack=JOBS bench=1
#                pack=JOBS caps the number of D4 or BigWig jobs running at
#                once. bench=1 keeps the timed scans one at a time.
# The CPUs really available stay at 256: the Slurm step the driver runs
# snakemake in is pinned to 256 cores. --cores and the thread counts only
# decide how many jobs Snakemake starts.

import sys
from pathlib import Path

# Paths in the config file are relative to the repo root (the parent of workflows/).
ROOT = Path(workflow.basedir).resolve().parents[1]
PATH_KEYS = [
    "manifest", "crams_dir", "reference", "chrom_sizes", "regions_bed",
    "pbz", "pbz_repo", "codecs", "results", "samples_file",
]


def from_root(path):
    path = Path(path).expanduser()
    return str(path if path.is_absolute() else ROOT / path)


CFG = {key: from_root(value) if key in PATH_KEYS else value for key, value in config.items()}
CORES = int(CFG["cores"])
RESULTS = CFG["results"]

# Every timed command must really get CORES threads.
if workflow.cores < CORES:
    raise ValueError(f"--cores must be at least {CORES}")

# The first n_samples manifest samples, in order, one per line.
with open(CFG["samples_file"]) as handle:
    SAMPLES = [line.strip() for line in handle if line.strip()]
if len(SAMPLES) != int(CFG["n_samples"]):
    raise ValueError(f"{CFG['samples_file']} has {len(SAMPLES)} samples, not n_samples={CFG['n_samples']}")


def cram(sample):
    return f"{CFG['crams_dir']}/{sample}.final.cram"


# The reference is a plain path, not a rule input: the driver copies it to local
# disk before the build, and every tool reads that copy.
REFERENCE = config["reference_local"]

# Where each tool's build leaves its data, inside the working directory.
PBZ_STORE = "data/cohort.pbz"
D4_FILE = "data/cohort.d4"
BIGWIGS = [f"data/{sample}.bw" for sample in SAMPLES]


def TIMED(command, timefile):
    """Wrap one command so /usr/bin/time writes wall s, user s, sys s and max RSS KB.

    `sync` is inside the timed command, so the time includes flushing the output
    to disk. Paths have no spaces, so they are not quoted inside sh -c '...'.
    """
    return f"/usr/bin/time -f '%e\\t%U\\t%S\\t%M' -o {timefile} sh -c '{command} && sync'"
