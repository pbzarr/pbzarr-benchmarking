#!/bin/bash
#SBATCH --job-name=part1
#SBATCH --partition=reserved
#SBATCH --account=vollger-reserved
#SBATCH --nodelist=rw244
#SBATCH --exclusive
#SBATCH --mem=1450000M
#SBATCH --time=1-00:00:00
#SBATCH --output=results/single/logs/part1-%j.log

# Part 1: time creating one depth track from one CRAM, for PBZ, D4 and BigWig.
# Submit from the repo root, after creating results/single/logs. To run the
# smoke test instead: CONFIG=config/single-smoke.yaml sbatch workflows/single/part1.sbatch
# If the job dies, submit it again: runs that already have a result file are skipped.
# The job asks for 1450000M, not --mem=0: all of rw244's memory (1532000M) is above
# the account's group limit (1487.5G), so Slurm rejects --mem=0 with AssocGrpMemLimit.
# --exclusive still gives the job the whole node.

set -euo pipefail

R=/scratch/ucgd/lustre-labs/vollger/users/cade/pbzarr-benchmarking
cd "$R"
export PATH="$R/.pixi/envs/default/bin:$PATH"
CONFIG=${CONFIG:-config/single.yaml}

# Read one value from the config file.
config_value() {
    python -c 'import sys, yaml; print(yaml.safe_load(open(sys.argv[1]))[sys.argv[2]])' "$CONFIG" "$1"
}
RESULTS=$(config_value results)
PBZ_REPO=$(config_value pbz_repo)
CODECS=$(config_value codecs)

# Record what ran, where and with which versions.
mkdir -p "$RESULTS"
{
    echo "hostname: $(hostname)"
    echo "slurm job id: $SLURM_JOB_ID"
    echo "date: $(date -Is)"
    echo "config: $CONFIG"
    echo "pbz commit: $(git -C "$PBZ_REPO" rev-parse --short HEAD)"
    echo "codecs: $(sha256sum "$CODECS")"
    echo "tool versions:"
    ls "$R/.pixi/envs/default/conda-meta" | grep -E '^(mosdepth|d4tools|megadepth|htslib)-'
    echo "squeue:"
    squeue -w "$(hostname)"
} > "$RESULTS/context.txt"

srun --ntasks=1 --cpus-per-task=256 --hint=nomultithread --mem=0 \
    snakemake -s workflows/single/Snakefile --configfile "$CONFIG" \
    --directory "${TMPDIR:?TMPDIR is not set}/single" --cores 256 \
    --resources bench=1 --keep-going
