#!/usr/bin/env bash
# Is d4tools merge slower when its files are on Lustre? Throwaway diagnostic.
#
# Makes per-sample D4 files with mosdepth for N tail-of-manifest samples, keeps
# one copy on Lustre and one on local disk, then times d4tools merge three ways
# with a cold node cache: Lustre in + out, Lustre in + local out, local in + out.
# Submit from the repo root:
#   sbatch workflows/checks/d4_merge_check.sh
#SBATCH --job-name=d4-merge-check
#SBATCH --partition=reserved
#SBATCH --account=vollger-reserved
#SBATCH --nodelist=rw244
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=results/checks/d4-merge-check-%j.log

set -euo pipefail

REPO=/scratch/ucgd/lustre-labs/vollger/users/cade/pbzarr-benchmarking
CRAMS=/scratch/ucgd/lustre-labs/vollger/data-sets/1kg-high-coverage/crams
REF=/scratch/ucgd/lustre-labs/vollger/assemblies/GRCh38_full_analysis_set_plus_decoy_hla.fa
export PATH=$REPO/.pixi/envs/default/bin:$PATH
EVICT="python $REPO/workflows/lib/evict.py"

N=${N:-32}
JOBS=6  # concurrent mosdepth jobs: about 4.5 GB each, inside the 32 GB limit
REPS=2

OUT=$REPO/results/checks/d4-merge/${SLURM_JOB_ID:-manual}
LUSTRE=$OUT/files
LOCAL=${TMPDIR:?run under Slurm, which sets TMPDIR to a per-job folder on /scratch/local}/d4-merge
mkdir -p "$LUSTRE" "$LOCAL/files"
trap 'rm -rf "$LUSTRE" "$LOCAL"' EXIT

mapfile -t SAMPLES < <(tail -n +2002 "$REPO/config/1kg-sr.tsv" | cut -f1 | head -"$N")
cp "$REF" "$REF.fai" "$LOCAL/"
ref_local=$LOCAL/$(basename "$REF")

# Per-sample D4 files, written locally and then copied to Lustre.
make_d4() {
    mosdepth -t 4 --fasta "$2" --d4 "$3/$1" "$4/$1.final.cram"
    mv "$3/$1.per-base.d4" "$3/$1.d4"
    rm -f "$3/$1".mosdepth.* "$3/$1".per-base.d4.csi
}
export -f make_d4
start=$(date +%s)
printf '%s\n' "${SAMPLES[@]}" | xargs -P "$JOBS" -I{} bash -c 'make_d4 "$@"' _ {} "$ref_local" "$LOCAL/files" "$CRAMS"
echo "mosdepth for $N samples: $(( $(date +%s) - start )) s"
cp "$LOCAL"/files/*.d4 "$LUSTRE/"
du -sh "$LUSTRE" "$LOCAL/files"

printf 'inputs\toutput\tn\trep\twall_s\tuser_s\tsys_s\tmax_rss_kb\tout_bytes\tstart\n' > "$OUT/merge.tsv"
merge() {  # merge INPUT_DIR OUTPUT_DIR LABEL_IN LABEL_OUT REP
    local in=$1 out=$2 sources=() s
    for s in "${SAMPLES[@]}"; do sources+=("$in/$s.d4:$s"); done
    rm -f "$out/merged.d4"
    $EVICT "$in"/*.d4
    local start; start=$(date -Is)
    /usr/bin/time -f '%e\t%U\t%S\t%M' -o "$LOCAL/time.tsv" d4tools merge "${sources[@]}" "$out/merged.d4"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$3" "$4" "$N" "$5" "$(cat "$LOCAL/time.tsv")" \
        "$(stat -c %s "$out/merged.d4")" "$start" >> "$OUT/merge.tsv"
    rm -f "$out/merged.d4"
}
for rep in $(seq "$REPS"); do
    merge "$LUSTRE" "$LUSTRE" lustre lustre "$rep"
    merge "$LUSTRE" "$LOCAL" lustre local "$rep"
    merge "$LOCAL/files" "$LOCAL" local local "$rep"
done

echo "done: $OUT"
