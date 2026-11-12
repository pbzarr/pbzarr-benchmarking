#!/bin/bash
# Submit part 2 for one config: one part2.sbatch job for each of its
# sample_counts. Every job asks for all of rw244, so Slurm runs them one at a time.
#   workflows/cohort/submit_part2.sh config/cohort-ont.yaml

set -euo pipefail
cd "$(dirname "$0")/../.."
config=${1:?usage: submit_part2.sh CONFIG}

mkdir -p results/logs
for n in $(.pixi/envs/default/bin/yq '.sample_counts[]' "$config"); do
    sbatch workflows/cohort/part2.sbatch "$config" "$n"
done
