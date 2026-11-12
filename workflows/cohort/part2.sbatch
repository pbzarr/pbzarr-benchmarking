#!/bin/bash
#SBATCH --job-name=part2
#SBATCH --partition=reserved
#SBATCH --account=vollger-reserved
#SBATCH --nodelist=rw244
#SBATCH --exclusive
#SBATCH --mem=1450000M
#SBATCH --time=3-00:00:00
#SBATCH --output=results/logs/part2-%j.log

# Part 2: build one cohort of the first N samples with each tool in turn, on the
# node's local disk, and time the scans on it. Submit from the repo root:
#   sbatch workflows/cohort/part2.sbatch config/cohort-ont.yaml 1019
# Results go to <results>/n<N>/. The job stops at once if that folder exists.
# (--mem=0 would ask for more memory than the account may use, so Slurm would
# reject it. --exclusive still gives the job the whole node.)

set -euo pipefail
CONFIG=${1:?usage: part2.sbatch CONFIG N}
N=${2:?usage: part2.sbatch CONFIG N}

R=/scratch/ucgd/lustre-labs/vollger/users/cade/pbzarr-benchmarking
cd "$R"
export PATH="$R/.pixi/envs/default/bin:$PATH"

CORES=$(yq .cores "$CONFIG")
REFERENCE=$(yq .reference "$CONFIG")
OUT=$(yq .results "$CONFIG")/n$N
WORK=$TMPDIR/part2

if [ -e "$OUT" ]; then
    echo "$OUT already exists. Delete it to run N=$N again." >&2
    exit 1
fi
mkdir -p "$OUT/timing" "$OUT/logs"
cp "$CONFIG" "$OUT/config.yaml"
{
    echo "start: $(date -Is) on $(hostname), Slurm job $SLURM_JOB_ID"
    echo "pbz: $($(yq .pbz "$CONFIG") --version)"
    echo "pbz commit: $(git -C "$(yq .pbz_repo "$CONFIG")" rev-parse --short HEAD)"
    ls .pixi/envs/default/conda-meta | grep -E '^(mosdepth|d4tools|megadepth|htslib)-'
} > "$OUT/context.txt"

# The scans start many threads and open many files, so raise both limits as far
# as they go.
ulimit -n "$(ulimit -Hn)"
ulimit -u "$(ulimit -Hu)"

# The samples: the first N of the manifest (its first column is sample_id). The
# config's `alignment` and `index` give each one's files, with {sample} in the path.
# The build is timed with all of these files cleared from the page cache.
cut -f1 "$(yq .manifest "$CONFIG")" | sed -n "2,$((N + 1))p" > "$OUT/samples.txt"
alignment=$(yq .alignment "$CONFIG")
index=$(yq .index "$CONFIG")
alignments=()
cold=()
while read -r sample; do
    alignments+=("${alignment//"{sample}"/$sample}")
    cold+=(--data "${alignment//"{sample}"/$sample}" --data "${index//"{sample}"/$sample}")
done < "$OUT/samples.txt"
read_bytes=$(du -cbL "${alignments[@]}" | tail -n 1 | cut -f1)

# CRAM input: copy the reference to local disk once, untimed. Every tool reads
# that copy. BAM input has `reference: null` and needs none.
reference=()
if [ "$REFERENCE" != null ]; then
    mkdir -p "$WORK"
    cp "$REFERENCE" "$WORK/reference.fa"
    cp "$REFERENCE.fai" "$WORK/reference.fa.fai"
    reference=("reference_local=$WORK/reference.fa")
fi

printf 'n\ttool\tthreads_per_job\tread_bytes\twall_s\tmax_rss_mib\tio_in\tcpu_s\tout_bytes\tout_files\n' > "$OUT/builds.tsv"

for tool in $(yq '.tools[]' "$CONFIG"); do
    # PBZ is one process on all cores. D4 and BigWig run one job per sample, all
    # at once, with the cores shared out: cores / N threads each, or 1 thread
    # each from N = cores up (then Snakemake runs `cores` jobs at a time).
    if [ "$tool" = pbz ]; then
        threads=$CORES
    else
        threads=$((N < CORES ? CORES / N : 1))
    fi
    args=(-s "workflows/cohort/$tool.smk" --cores "$CORES" --directory "$WORK/$tool"
        --configfile "$CONFIG" --config "n_samples=$N" "samples_file=$OUT/samples.txt"
        "results=$OUT" "threads_per_job=$threads" "${reference[@]}")

    # The build: one timed run. The time includes the `sync` that flushes the
    # output to disk.
    python workflows/create/timed_trials.py --cache cold --trials 1 \
        --output "$OUT/timing/build-$tool.tsv" "${cold[@]}" \
        -- sh -c 'snakemake build "$@" && sync' sh "${args[@]}"

    # The build's row. Columns 3-6 of the timing file are wall s, peak RSS MiB,
    # MiB read from storage and CPU s. Snakemake leaves an empty
    # .snakemake_timestamp file in folder outputs, which is not counted.
    timing=$(tail -n 1 "$OUT/timing/build-$tool.tsv" | cut -f3-6)
    out_bytes=$(du -sb --exclude=.snakemake_timestamp "$WORK/$tool/data" | cut -f1)
    out_files=$(find "$WORK/$tool/data" -type f ! -name .snakemake_timestamp | wc -l)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$N" "$tool" "$threads" "$read_bytes" \
        "$timing" "$out_bytes" "$out_files" >> "$OUT/builds.tsv"

    # The scans, one at a time, each with all cores. If they fail or run longer
    # than 16 hours, failed.txt says so and the job goes on to the next tool.
    timeout --kill-after=60 16h snakemake scans --resources bench=1 "${args[@]}" \
        || echo "$tool scans: exit status $? (124 means timed out)" >> "$OUT/failed.txt"

    if [ -d "$WORK/$tool/logs" ]; then
        cp -r "$WORK/$tool/logs" "$OUT/logs/$tool"
    fi
    rm -rf "${WORK:?}/$tool"
done

date -Is > "$OUT/DONE"
