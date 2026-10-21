#!/bin/bash
#SBATCH --job-name=part2
#SBATCH --partition=reserved
#SBATCH --account=vollger-reserved
#SBATCH --nodelist=rw244
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=1-00:00:00
#SBATCH --output=results/cohort/logs/part2-%j.log

# Part 2: build one cohort of N samples with PBZ, D4 and BigWig, one tool at a
# time. Usage: part2.sbatch N BUILD   (submit_part2.sh submits every job)
# Each build is its own Slurm step, named build-<tool>, so sacct gives its wall
# time, CPU time and max RSS. The scans run on build 1 only. Then the tool's
# data is deleted. Submit with CONFIG= and PACKING= set to use other files, for
# example the smoke test: config/cohort-smoke.yaml and config/packing-smoke.yaml.
# If the job dies, submit it again: a job without $OUT/DONE starts over, one
# with it exits.

set -euo pipefail

N=${1:?usage: part2.sbatch N BUILD}
BUILD=${2:?usage: part2.sbatch N BUILD}

R=/scratch/ucgd/lustre-labs/vollger/users/cade/pbzarr-benchmarking
cd "$R"
export PATH="$R/.pixi/envs/default/bin:$PATH"
CONFIG=${CONFIG:-config/cohort.yaml}
PACKING=${PACKING:-config/packing.yaml}

# Read one value from the config file, one line per list item.
config_value() {
    python -c '
import sys, yaml
value = yaml.safe_load(open(sys.argv[1]))[sys.argv[2]]
for item in value if isinstance(value, list) else [value]:
    print(item)
' "$CONFIG" "$1"
}
RESULTS=$(config_value results)
case "$RESULTS" in /*) ;; *) RESULTS=$R/$RESULTS ;; esac
MANIFEST=$(config_value manifest)
CRAMS_DIR=$(config_value crams_dir)
REFERENCE=$(config_value reference)
CODECS=$(config_value codecs)
PBZ_REPO=$(config_value pbz_repo)
CORES=$(config_value cores)
DISK_GB=$(config_value disk_per_sample_gb)

OUT=$RESULTS/n$N/b$BUILD
WORK=${TMPDIR:?TMPDIR is not set}/cohort

# Already done: nothing to do.
if [ -e "$OUT/DONE" ]; then
    echo "already done"
    exit 0
fi
[ -f "$PACKING" ] || { echo "no packing file: $PACKING" >&2; exit 1; }
rm -rf "$OUT"
mkdir -p "$OUT"

# Keep each tool's logs on Lustre, also when the job fails. Called after every
# tool and from the EXIT trap.
save_logs() {
    local dir tool
    for dir in "$WORK"/*/logs; do
        [ -d "$dir" ] || continue
        tool=$(basename "$(dirname "$dir")")
        mkdir -p "$OUT/logs/$tool"
        cp -r "$dir/." "$OUT/logs/$tool/" || true
    done
}
trap save_logs EXIT

# Record what ran, where and with which versions.
{
    echo "hostname: $(hostname)"
    echo "slurm job id: $SLURM_JOB_ID"
    echo "date: $(date -Is)"
    echo "config: $CONFIG"
    echo "packing: $PACKING"
    echo "n: $N"
    echo "build: $BUILD"
    echo "pbz commit: $(git -C "$PBZ_REPO" rev-parse --short HEAD || true)"
    echo "codecs: $(sha256sum "$CODECS")"
    echo "tool versions:"
    ls "$R/.pixi/envs/default/conda-meta" | grep -E '^(mosdepth|d4tools|megadepth|htslib)-' || true
    echo "squeue:"
    squeue -w "$(hostname)" || true
} > "$OUT/context.txt"

# Open files: raise the limit as far as it goes, and stop if it is below 4 * N + 256.
ulimit -n "$(ulimit -Hn)" || true
limit=$(ulimit -n)
if [ "$limit" != unlimited ] && [ "$limit" -lt $((4 * N + 256)) ]; then
    echo "open-file limit $limit is below $((4 * N + 256))" >&2
    exit 1
fi

# Disk: the node-local space under TMPDIR must be larger than N * disk_per_sample_gb GB.
free_bytes=$(df -Pk "$TMPDIR" | awk 'NR == 2 {printf "%.0f\n", $4 * 1024}')
need_bytes=$(awk -v n="$N" -v gb="$DISK_GB" 'BEGIN {printf "%.0f\n", n * gb * 1e9}')
if [ "$free_bytes" -le "$need_bytes" ]; then
    echo "not enough local disk: $free_bytes bytes free under $TMPDIR, need more than $need_bytes" >&2
    exit 1
fi

# Setup. The reference is copied once, untimed; every tool reads the local copy.
mkdir -p "$WORK/ref"
cp "$REFERENCE" "$WORK/ref/reference.fa"
cp "$REFERENCE.fai" "$WORK/ref/reference.fa.fai"

# The first N samples of the manifest.
awk -F'\t' -v n="$N" '
    NR == 1 { for (i = 1; i <= NF; i++) if ($i == "sample_id") col = i; if (!col) exit 1; next }
    NR > n + 1 { exit }
    { print $col }
' "$MANIFEST" > "$OUT/samples.txt"
if [ "$(wc -l < "$OUT/samples.txt")" -ne "$N" ]; then
    echo "$MANIFEST has fewer than $N samples" >&2
    exit 1
fi

# Bytes the tools read: the CRAMs of these samples.
READ_BYTES=$(python -c '
import os, sys
crams_dir, samples = sys.argv[1:3]
print(sum(os.path.getsize(f"{crams_dir}/{line.strip()}.final.cram") for line in open(samples) if line.strip()))
' "$CRAMS_DIR" "$OUT/samples.txt")

# Every CRAM and .crai, to clear from the page cache before each build.
awk -v dir="$CRAMS_DIR" '{ print dir "/" $1 ".final.cram"; print dir "/" $1 ".final.cram.crai" }' \
    "$OUT/samples.txt" > "$WORK/evict.txt"

# Tool order for this job: the config's tools, shuffled with the fixed seed
# "$N-$BUILD". The seed is hashed first: shuf reads only the first few bytes of
# its random source, so "256-1", "256-2" and "256-3" alone give the same order.
seed=$(printf '%s' "$N-$BUILD" | sha256sum | cut -d' ' -f1)
config_value tools | shuf --random-source=<(yes "$seed") > "$OUT/order.txt"
ORDER=$(cat "$OUT/order.txt")

table=$OUT/builds.tsv
position=0
for tool in $ORDER; do
    position=$((position + 1))
    data=$WORK/$tool/data

    # Cold inputs for the build.
    xargs python workflows/lib/evict.py < "$WORK/evict.txt"

    # Threads per job and jobs at once. PBZ is one process using every core.
    case "$tool" in
        pbz) THREADS=$CORES; JOBS=1 ;;
        d4 | bigwig)
            packing=$(python workflows/cohort/packing.py for "$PACKING" "$tool" "$N")
            read -r THREADS JOBS <<< "$packing"
            ;;
        *) echo "unknown tool: $tool" >&2; exit 1 ;;
    esac

    # Snakemake starts jobs by their declared threads, so --cores must cover
    # JOBS * THREADS. The Slurm steps below keep the real CPUs at 256.
    snakemake_cores=$((JOBS * THREADS > CORES ? JOBS * THREADS : CORES))
    snakemake_config=(--directory "$WORK/$tool" --configfile "$CONFIG"
        --config "n_samples=$N" "samples_file=$OUT/samples.txt" "results=$OUT" "build=$BUILD"
        "threads_per_job=$THREADS" "reference_local=$WORK/ref/reference.fa")

    # The timed build: its own Slurm step, so sacct reports it.
    srun --ntasks=1 --cpus-per-task=256 --hint=nomultithread --mem=0 \
        --acctg-freq=task=10 --job-name="build-$tool" \
        sh -c 'snakemake "$@" && sync' sh \
        -s "workflows/cohort/$tool.smk" build --cores "$snakemake_cores" \
        --resources "pack=$JOBS" bench=1 "${snakemake_config[@]}"

    # The step's accounting row, once it is COMPLETED (up to 60 s).
    row=
    for _ in $(seq 30); do
        row=$(sacct -j "$SLURM_JOB_ID" -n -P --units=K \
            -o JobID,JobName,State,ElapsedRaw,TotalCPU,UserCPU,SystemCPU,MaxRSS \
            | awk -F'|' -v name="build-$tool" '$2 == name && $3 == "COMPLETED"' || true)
        [ -n "$row" ] && break
        sleep 2
    done
    [ -n "$row" ] || { echo "no COMPLETED sacct row for build-$tool" >&2; exit 1; }
    IFS='|' read -r _ _ _ elapsed_s total_cpu user_cpu system_cpu max_rss <<< "$row"

    # Output size, without Snakemake's .snakemake_timestamp file. du -sb counts
    # directory entries too, so a build with no non-empty file is also refused.
    out_bytes=$(du -sb --exclude=.snakemake_timestamp "$data" | cut -f1)
    out_files=$(find "$data" -type f ! -name .snakemake_timestamp | awk 'END {print NR}')
    nonempty_files=$(find "$data" -type f -size +0c ! -name .snakemake_timestamp | awk 'END {print NR}')
    if [ "$out_bytes" -eq 0 ] || [ "$nonempty_files" -eq 0 ]; then
        echo "$tool build wrote no output bytes" >&2
        exit 1
    fi

    # Most disk the build held: D4 keeps every per-sample file next to the merged one.
    peak_disk_bytes=$out_bytes
    if [ "$tool" = d4 ]; then
        peak_disk_bytes=$((out_bytes + $(< "$OUT/d4_per_sample_bytes.txt")))
    fi

    if [ ! -e "$table" ]; then
        printf 'n\tbuild\ttool\torder\tthreads_per_job\tjobs\tread_bytes\telapsed_s\ttotal_cpu\tuser_cpu\tsystem_cpu\tmax_rss_kb\tout_bytes\tout_files\tpeak_disk_bytes\n' > "$table"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$N" "$BUILD" "$tool" "$position" "$THREADS" "$JOBS" "$READ_BYTES" "$elapsed_s" \
        "$total_cpu" "$user_cpu" "$system_cpu" "${max_rss%K}" "$out_bytes" "$out_files" \
        "$peak_disk_bytes" >> "$table"

    # The scans, on build 1 only. Same pinning and config as the build, so
    # Snakemake sees nothing to rebuild.
    if [ "$BUILD" = 1 ]; then
        srun --ntasks=1 --cpus-per-task=256 --hint=nomultithread --mem=0 \
            --job-name="scans-$tool" \
            snakemake -s "workflows/cohort/$tool.smk" scans --cores "$snakemake_cores" \
            --resources bench=1 "${snakemake_config[@]}"
    fi

    save_logs
    rm -rf "${WORK:?}/$tool"
done

touch "$OUT/DONE"
