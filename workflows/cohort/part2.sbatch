#!/bin/bash
#SBATCH --job-name=part2
#SBATCH --partition=reserved
#SBATCH --account=vollger-reserved
#SBATCH --nodelist=rw244
#SBATCH --exclusive
#SBATCH --mem=1450000M
#SBATCH --time=1-00:00:00
#SBATCH --output=results/cohort/logs/part2-%j.log

# Part 2: build one cohort of N samples with PBZ, D4 and BigWig, one tool at a
# time. Usage: part2.sbatch N BUILD   (submit_part2.sh submits every job)
# Each build is timed by workflows/create/timed_trials.py, the tool the scans use
# too: one trial with a cold page cache (it clears every CRAM and .crai first).
# It writes $OUT/timing/build-<tool>.tsv, and builds.tsv gets one row per tool
# from it. The scans run on build 1 only. Then the tool's data is deleted.
# Submit with CONFIG= set to use another file, for example the smoke test:
#   CONFIG=config/cohort-smoke.yaml sbatch workflows/cohort/part2.sbatch 4 1
# If the job dies, submit it again: a job without $OUT/DONE starts over, one
# with it exits.
# The job asks for 1450000M, not --mem=0: all of rw244's memory (1532000M) is above
# the account's group limit (1487.5G), so Slurm rejects --mem=0 with AssocGrpMemLimit.
# --exclusive still gives the job the whole node.
#
# Threads per job (THREADS) and jobs at once (JOBS) for each build:
#   PBZ                   THREADS = cores, JOBS = 1.
#   D4, BigWig, N < cores JOBS = N, THREADS = max(1, floor(cores / N)).
#   D4, BigWig, N >= cores THREADS = packing.<tool>.threads and JOBS =
#                         min(packing.<tool>.jobs, N), from the config file, set
#                         by hand from part 1's results. The job refuses to
#                         start if they are null. builds.tsv records the JOBS used.
# The build runs snakemake with --cores JOBS * THREADS, the scans with
# --cores <cores>. See workflows/cohort/common.smk.

set -euo pipefail

N=${1:?usage: part2.sbatch N BUILD}
BUILD=${2:?usage: part2.sbatch N BUILD}

R=/scratch/ucgd/lustre-labs/vollger/users/cade/pbzarr-benchmarking
cd "$R"
export PATH="$R/.pixi/envs/default/bin:$PATH"
CONFIG=${CONFIG:-config/cohort.yaml}

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

# Set THREADS (threads per job) and JOBS (jobs at once) for one tool. Stops the
# whole job if a value that is needed is missing from the config.
pack() {
    local tool=$1 threads jobs
    case "$tool" in
        pbz)
            # One process using every core.
            THREADS=$CORES
            JOBS=1
            ;;
        d4 | bigwig)
            if [ "$N" -lt "$CORES" ]; then
                # Fewer samples than cores: every sample at once, the cores shared out evenly.
                JOBS=$N
                THREADS=$((CORES / N > 1 ? CORES / N : 1))
            else
                threads=$(config_value "packing.$tool.threads" || true)
                jobs=$(config_value "packing.$tool.jobs" || true)
                if ! [[ $threads =~ ^[1-9][0-9]*$ && $jobs =~ ^[1-9][0-9]*$ ]]; then
                    echo "N=$N is not below cores=$CORES, so $CONFIG must give packing.$tool.threads and packing.$tool.jobs as whole numbers (now: threads='$threads', jobs='$jobs'). Set them from part 1's results." >&2
                    exit 1
                fi
                THREADS=$threads
                # Never more jobs than samples.
                JOBS=$((jobs < N ? jobs : N))
            fi
            ;;
        *) echo "unknown tool: $tool" >&2; exit 1 ;;
    esac
}
# Check every tool in the config now, before anything is deleted or built, so a
# missing packing value stops the job at once.
for tool in $(config_value tools); do
    pack "$tool"
done

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
    echo "n: $N"
    echo "build: $BUILD"
    echo "pbz commit: $(git -C "$PBZ_REPO" rev-parse --short HEAD || true)"
    echo "codecs: $(sha256sum "$CODECS")"
    echo "tool versions:"
    ls "$R/.pixi/envs/default/conda-meta" | grep -E '^(mosdepth|d4tools|megadepth|htslib)-' || true
    echo "squeue:"
    squeue -w "$(hostname)" || true
    echo "lscpu:"
    lscpu || true
} > "$OUT/context.txt"

# Open files: raise the limit as far as it goes, and stop if it is below 4 * N + 256.
ulimit -n "$(ulimit -Hn)" || true
# Threads: the scans start many threads (a pool of 256 worker processes, each with
# threads of its own). Slurm's default limit is 16384 threads per user, and a
# scan that reaches it hangs. Raise it as far as it goes.
ulimit -u "$(ulimit -Hu)" || true
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

# Every CRAM and .crai as a --data argument of timed_trials.py, which clears them
# from the page cache before the build.
data_args=()
while read -r sample; do
    data_args+=(--data "$CRAMS_DIR/$sample.final.cram" --data "$CRAMS_DIR/$sample.final.cram.crai")
done < "$OUT/samples.txt"

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
    pack "$tool"

    # Snakemake starts jobs by their declared threads, so --cores = JOBS * THREADS
    # runs exactly JOBS per-sample jobs at once (for PBZ it is just the cores).
    snakemake_cores=$((JOBS * THREADS))
    snakemake_config=(--directory "$WORK/$tool" --configfile "$CONFIG"
        --config "n_samples=$N" "samples_file=$OUT/samples.txt" "results=$OUT" "build=$BUILD"
        "threads_per_job=$THREADS" "reference_local=$WORK/ref/reference.fa")

    # The timed build: one trial, CRAMs and .crai cold. `sync` is inside the
    # timed command, so the time includes flushing the output to disk.
    timing=$OUT/timing/build-$tool.tsv
    python workflows/create/timed_trials.py --cache cold --trials 1 --output "$timing" \
        "${data_args[@]}" \
        -- sh -c 'snakemake "$@" && sync' sh \
        -s "workflows/cohort/$tool.smk" build --cores "$snakemake_cores" "${snakemake_config[@]}"

    # The build's row in the timing table: wall s, CPU s (user + sys), peak RSS in
    # MiB (summed over the process tree, sampled) and MiB read from storage.
    row=$(awk -F'\t' '
        NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        NR == 2 { print $col["s"] "\t" $col["cpu_time"] "\t" $col["max_rss"] "\t" $col["io_in"] }
    ' "$timing")
    [ -n "$row" ] || { echo "no timing row in $timing" >&2; exit 1; }
    IFS=$'\t' read -r wall_s cpu_s max_rss_mib io_in <<< "$row"

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
        printf 'n\tbuild\ttool\torder\tthreads_per_job\tjobs\tread_bytes\twall_s\tcpu_s\tmax_rss_mib\tio_in\tout_bytes\tout_files\tpeak_disk_bytes\n' > "$table"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$N" "$BUILD" "$tool" "$position" "$THREADS" "$JOBS" "$READ_BYTES" "$wall_s" \
        "$cpu_s" "$max_rss_mib" "$io_in" "$out_bytes" "$out_files" "$peak_disk_bytes" >> "$table"

    # The scans, on build 1 only. Same config as the build, so Snakemake sees
    # nothing to rebuild. They run one at a time (bench=1), each with the cores.
    if [ "$BUILD" = 1 ]; then
        snakemake -s "workflows/cohort/$tool.smk" scans --cores "$CORES" \
            --resources bench=1 "${snakemake_config[@]}"
    fi

    save_logs
    rm -rf "${WORK:?}/$tool"
done

touch "$OUT/DONE"
