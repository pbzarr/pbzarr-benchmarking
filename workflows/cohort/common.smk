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
# Snakemake schedules by the threads each rule declares, so the driver passes
# --cores like this:
#   D4 and BigWig build   --cores JOBS * THREADS. Each per-sample job declares
#                         THREADS (threads_per_job), so Snakemake runs exactly
#                         JOBS of them at once.
#   PBZ build, scans      --cores <cores> (the config's `cores`).
#   scans                 --resources bench=1, so the timed scans run one at a
#                         time.
# The PBZ import and every scan declare threads = cores, and their shell starts
# with `test {threads} -eq <cores>`. If Snakemake ever lowered their threads
# (because --cores was too small), the rule fails instead of running a timed
# command with fewer threads.
#
# The build is timed from outside by the driver, with timed_trials.py. The scans
# time themselves with the same tool, inside the rules. The one step timed
# inside a workflow is D4's merge (see d4.smk).

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


# What the scans in all three workflows use. timed_trials.py is the same timing
# tool the creation benchmark in workflows/create uses; the scan scripts are in
# the same folder and are used in place.
CREATE = ROOT / "workflows/create"
TIMED_TRIALS = str(CREATE / "timed_trials.py")
PYTHON = sys.executable
REGIONS_BED = CFG["regions_bed"]
TRIALS = int(CFG["scan_trials"])

wildcard_constraints:
    analysis="full_scan|region_scan"


# The whole-genome region list for the full scans: every contig, 0 to its length.
rule full_scan_bed:
    input:
        CFG["chrom_sizes"]
    output:
        "full_scan.bed"
    shell:
        r"""awk 'BEGIN {{OFS="\t"}} $2 > 0 {{print $1, 0, $2}}' {input:q} > {output:q}"""
