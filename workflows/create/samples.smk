import csv
from pathlib import Path


WORKFLOW_DIR = Path(workflow.basedir).resolve()
ROOT = WORKFLOW_DIR.parents[1]
MANIFEST = str(ROOT / config["manifest"])
REFERENCE = str(ROOT / config["reference"])
REFERENCE_INDEX = str(ROOT / config.get("reference_index", config["reference"] + ".fai"))
CHROM_SIZES = str(ROOT / config["chrom_sizes"])

with open(MANIFEST) as handle:
    SAMPLES = list(csv.DictReader(handle, delimiter="\t"))
for sample in SAMPLES:
    for field in ("alignment", "index"):
        sample[field] = str(ROOT / sample[field])


def selected_samples(n):
    n = int(n)
    if n < 1 or n > len(SAMPLES):
        raise ValueError(f"Requested {n} samples, but {MANIFEST} has {len(SAMPLES)}")
    return SAMPLES[:n]
