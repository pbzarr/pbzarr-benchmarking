#!/bin/bash
# Submit part 2: one job per cohort size (build 1), then the noise-check builds
# 2 ... noise_check.builds at noise_check.n. The jobs run one after another: each
# starts when the one before it has ended, however it ended. Run from anywhere.
# To use another config file, for example the smoke test:
#   CONFIG=config/cohort-smoke.yaml workflows/cohort/submit_part2.sh

set -euo pipefail

R=$(cd "$(dirname "$0")/../.." && pwd)
cd "$R"
export PATH="$R/.pixi/envs/default/bin:$PATH"
# part2.sbatch reads this from its environment, which sbatch copies from here.
export CONFIG=${CONFIG:-config/cohort.yaml}

# Read one value from the config file: a dotted name goes down into mappings
# (packing.d4.jobs), a list gives one line per item, and null gives an empty line.
config_value() {
    python -c '
import sys, yaml
value = yaml.safe_load(open(sys.argv[1]))
for key in sys.argv[2].split("."):
    value = value[key]
for item in value if isinstance(value, list) else [value]:
    print("" if item is None else item)
' "$CONFIG" "$1"
}
RESULTS=$(config_value results)

# Refuse now, before anything is submitted, if a job would need a packing value
# that is not set. part2.sbatch needs packing.<tool>.threads and .jobs for D4 and
# BigWig whenever a sample count is not below cores, and it stops at once without
# them, which would waste the jobs that ran before it.
cores=$(config_value cores)
noise_n=$(config_value noise_check.n)
noise_builds=$(config_value noise_check.builds)
sizes="$(config_value sample_counts) $noise_n"
for tool in $(config_value tools); do
    case "$tool" in d4 | bigwig) ;; *) continue ;; esac
    threads=$(config_value "packing.$tool.threads" || true)
    jobs=$(config_value "packing.$tool.jobs" || true)
    [[ $threads =~ ^[1-9][0-9]*$ && $jobs =~ ^[1-9][0-9]*$ ]] && continue
    for n in $sizes; do
        if [ "$n" -ge "$cores" ]; then
            echo "Nothing submitted. N=$n is not below cores=$cores, so $CONFIG must give packing.$tool.threads and packing.$tool.jobs as whole numbers (now: threads='$threads', jobs='$jobs'). Set them from part 1's results." >&2
            exit 1
        fi
    done
done

mkdir -p "$RESULTS/logs"

# Submit part2.sbatch N BUILD, after the job submitted before it. Sets $previous.
previous=
submit() {
    local n=$1 build=$2
    local args=(--parsable --output="$RESULTS/logs/part2-%j.log")
    if [ -n "$previous" ]; then
        args+=(--dependency="afterany:$previous")
    fi
    previous=$(sbatch "${args[@]}" workflows/cohort/part2.sbatch "$n" "$build")
    previous=${previous%%;*}
    echo "n=$n build=$build job=$previous"
}

for n in $(config_value sample_counts); do
    submit "$n" 1
done
for ((build = 2; build <= noise_builds; build++)); do
    submit "$noise_n" "$build"
done
