#!/usr/bin/env bash
# Lustre and page-cache checks for the creation benchmark. Throwaway diagnostics.
#
#   1. cache:  can posix_fadvise drop a Lustre CRAM from this node's page cache?
#   2. reads:  aggregate read speed from Lustre vs node-local disk, 1-64 streams.
#   3. tools:  pbz, mosdepth and Megadepth on one sample, CRAM on Lustre vs local.
#
# Uses CRAMs from the tail of the manifest, so the first 1024 benchmark samples
# stay untouched. Submit from the repo root:
#   mkdir -p results/checks && sbatch workflows/checks/lustre_check.sh
#SBATCH --job-name=lustre-check
#SBATCH --partition=reserved
#SBATCH --account=vollger-reserved
#SBATCH --nodelist=rw244
#SBATCH --cpus-per-task=64
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=results/checks/lustre-check-%j.log

set -euo pipefail

REPO=/scratch/ucgd/lustre-labs/vollger/users/cade/pbzarr-benchmarking
CRAMS=/scratch/ucgd/lustre-labs/vollger/data-sets/1kg-high-coverage/crams
REF=/scratch/ucgd/lustre-labs/vollger/assemblies/GRCh38_full_analysis_set_plus_decoy_hla.fa
PBZ=$HOME/pbzarr/target/release/pbz
export PATH=$REPO/.pixi/envs/default/bin:$PATH
EVICT="python $REPO/workflows/lib/evict.py"

OUT=$REPO/results/checks/lustre/${SLURM_JOB_ID:-manual}
LOCAL=${TMPDIR:?run under Slurm, which sets TMPDIR to a per-job folder on /scratch/local}/lustre-check
mkdir -p "$OUT" "$LOCAL/reads" "$LOCAL/out"
trap 'rm -rf "$LOCAL"' EXIT

READ_BLOCKS=1024  # 4 MiB blocks: read the first 4 GiB of each CRAM
REPS=2

# Samples 2001+ of the manifest: 1 for the cache check, 105 for reads, 1 for tools.
mapfile -t SAMPLES < <(tail -n +2002 "$REPO/config/1kg-sr.tsv" | cut -f1 | head -107)
cram() { echo "$CRAMS/$1.final.cram"; }
resident() { fincore --bytes --noheadings --output RES "$@" | awk '{s += $1} END {print s + 0}'; }
now() { date +%s.%N; }

{
    hostname; nproc; lscpu | grep -E 'Model name|Socket|Core|Thread|NUMA node\(s\)'
    df -h /scratch/local "$CRAMS"
    lfs getstripe -c -S "$(cram "${SAMPLES[0]}")"
    lctl get_param llite.*.max_read_ahead_mb llite.*.max_read_ahead_per_file_mb
    "$PBZ" --version; mosdepth --version; megadepth --version 2>&1 | head -1
} > "$OUT/node.txt" 2>&1 || true

# Copy the tool-check sample and the reference to local disk first, so the Lustre
# servers have done other work before the timed Lustre runs read this CRAM again.
sample=${SAMPLES[106]}
cp "$(cram "$sample")" "$(cram "$sample").crai" "$REF" "$REF.fai" "$LOCAL/"

# 1. cache: evict, read cold, read warm, evict again, read cold again.
f=$(cram "${SAMPLES[0]}")
size=$(stat -c %s "$f")
printf 'step\tresident_bytes\tfile_bytes\twall_s\tgib_per_s\n' > "$OUT/cache.tsv"
timed_cat() {
    local start end
    start=$(now); cat "$f" > /dev/null; end=$(now)
    awk -v s="$start" -v e="$end" -v b="$size" 'BEGIN {printf "%.2f\t%.3f", e - s, b / 2^30 / (e - s)}'
}
$EVICT "$f"
printf 'after_evict\t%s\t%s\t\t\n' "$(resident "$f")" "$size" >> "$OUT/cache.tsv"
t=$(timed_cat); printf 'cold_read\t%s\t%s\t%s\n' "$(resident "$f")" "$size" "$t" >> "$OUT/cache.tsv"
t=$(timed_cat); printf 'warm_read\t%s\t%s\t%s\n' "$(resident "$f")" "$size" "$t" >> "$OUT/cache.tsv"
$EVICT "$f"
printf 'after_evict_2\t%s\t%s\t\t\n' "$(resident "$f")" "$size" >> "$OUT/cache.tsv"
t=$(timed_cat); printf 'cold_read_2\t%s\t%s\t%s\n' "$(resident "$f")" "$size" "$t" >> "$OUT/cache.tsv"

# 2. reads: P parallel streams, each reading the first 4 GiB of its own file.
printf 'source\tstreams\tgib\twall_s\tgib_per_s\n' > "$OUT/reads.tsv"
parallel_read() {  # parallel_read SOURCE FILE...
    local source=$1; shift
    $EVICT "$@"
    local start end
    start=$(now)
    for file in "$@"; do dd if="$file" of=/dev/null bs=4M count=$READ_BLOCKS status=none & done
    wait
    end=$(now)
    awk -v src="$source" -v n=$# -v s="$start" -v e="$end" -v blocks=$READ_BLOCKS \
        'BEGIN {gib = n * blocks * 4 / 1024; printf "%s\t%d\t%d\t%.2f\t%.3f\n", src, n, gib, e - s, gib / (e - s)}' \
        >> "$OUT/reads.tsv"
}
next=1  # distinct Lustre files per level, so the Lustre servers' own cache can't help
lustre64=()
for streams in 1 8 32 64; do
    files=()
    for ((i = 0; i < streams; i++)); do files+=("$(cram "${SAMPLES[next + i]}")"); done
    next=$((next + streams))
    parallel_read lustre "${files[@]}"
    lustre64=("${files[@]}")
done

# Copy the first 4 GiB of the last 64 files to local disk, timed as one parallel copy.
start=$(now)
for file in "${lustre64[@]}"; do
    dd if="$file" of="$LOCAL/reads/$(basename "$file")" bs=4M count=$READ_BLOCKS status=none &
done
wait
end=$(now)
awk -v s="$start" -v e="$end" -v blocks=$READ_BLOCKS \
    'BEGIN {gib = 64 * blocks * 4 / 1024; printf "copy_to_local\t64\t%d\t%.2f\t%.3f\n", gib, e - s, gib / (e - s)}' \
    >> "$OUT/reads.tsv"
mapfile -t locals < <(ls "$LOCAL"/reads/*)
for streams in 1 8 32 64; do parallel_read local "${locals[@]:0:streams}"; done
rm -rf "$LOCAL/reads"

# 3. tools: one sample, CRAM on Lustre vs a local copy. The reference is local for both.
ref_local=$LOCAL/$(basename "$REF")
declare -A SRC=([lustre]=$(cram "$sample") [local]=$LOCAL/$(basename "$(cram "$sample")"))

printf 'tool\tthreads\tsource\trep\twall_s\tuser_s\tsys_s\tmax_rss_kb\tstart\n' > "$OUT/tools.tsv"
run_tool() {  # run_tool TOOL THREADS SOURCE REP
    local tool=$1 threads=$2 source=$3 rep=$4 input=${SRC[$3]} out=$LOCAL/out/run
    rm -rf "$out"; mkdir -p "$out"
    $EVICT "$input" "$input.crai"
    local start; start=$(date -Is)
    case $tool in
        pbz) set -- "$PBZ" import bam "$input:$sample" --reference "$ref_local" --threads "$threads" \
                 --track depth --no-progress -o "$out/cohort.pbz" ;;
        mosdepth) set -- mosdepth -t "$threads" --fasta "$ref_local" --d4 "$out/coverage" "$input" ;;
        megadepth) set -- megadepth "$input" --fasta "$ref_local" --bigwig --filter-out 1796 \
                       --threads "$threads" --prefix "$out/sample" ;;
    esac
    if ! /usr/bin/time -f '%e\t%U\t%S\t%M' -o "$out/time.tsv" "$@" > "$out/stdout.log" 2> "$out/stderr.log"; then
        cp "$out/stderr.log" "$OUT/failed-$tool-$threads-$source-$rep.log"
        printf '%s\t%s\t%s\t%s\tfailed\t\t\t\t%s\n' "$tool" "$threads" "$source" "$rep" "$start" >> "$OUT/tools.tsv"
        return 0
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$tool" "$threads" "$source" "$rep" "$(cat "$out/time.tsv")" "$start" \
        >> "$OUT/tools.tsv"
}
for rep in $(seq "$REPS"); do
    for spec in "pbz 64" "pbz 8" "mosdepth 4" "megadepth 4"; do
        for source in lustre local; do
            run_tool $spec "$source" "$rep"
        done
    done
done

echo "done: $OUT"
