#!/bin/bash
#SBATCH --job-name=part1
#SBATCH --partition=reserved
#SBATCH --account=vollger-reserved
#SBATCH --nodelist=rw244
#SBATCH --exclusive
#SBATCH --mem=1450000M
#SBATCH --time=2-00:00:00
#SBATCH --output=results/logs/part1-%j.log

# Part 1: time creating one depth track from one alignment file with each tool,
# then check where the tools agree. One config per data type. Submit from the
# repo root:
#   sbatch workflows/single/part1.sbatch config/single-ont.yaml
# If the job dies, submit it again: runs that already have a result file are skipped.
# (--mem=0 would ask for more memory than the account may use, so Slurm would
# reject it. --exclusive still gives the job the whole node.)

set -euo pipefail
CONFIG=${1:?usage: part1.sbatch CONFIG}

R=/scratch/ucgd/lustre-labs/vollger/users/cade/pbzarr-benchmarking
cd "$R"
export PATH="$R/.pixi/envs/default/bin:$PATH"

RESULTS=$(yq .results "$CONFIG")
mkdir -p "$RESULTS"
# A copy of the config. If it changed since the last job, the old copy is kept
# as config.yaml.~1~.
cmp -s "$CONFIG" "$RESULTS/config.yaml" || cp --backup=numbered "$CONFIG" "$RESULTS/config.yaml"
{
    echo "start: $(date -Is) on $(hostname), Slurm job $SLURM_JOB_ID"
    echo "pbz: $($(yq .pbz "$CONFIG") --version)"
    echo "pbz commit: $(git -C "$(yq .pbz_repo "$CONFIG")" rev-parse --short HEAD)"
    ls .pixi/envs/default/conda-meta | grep -E '^(mosdepth|d4tools|megadepth|htslib)-'
} >> "$RESULTS/context.txt"

snakemake -s workflows/single/Snakefile --configfile "$CONFIG" \
    --directory "$TMPDIR/single" --cores 256 --resources bench=1 --keep-going
