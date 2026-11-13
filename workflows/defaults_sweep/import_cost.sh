#!/bin/bash
#SBATCH --job-name=pbz-import-cost
#SBATCH --partition=reserved
#SBATCH --account=vollger-reserved
#SBATCH -n 1
#SBATCH -c 224
#SBATCH --mem=880G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/ucgd/lustre-labs/vollger/users/cade/defaults-sweep/logs/import-cost-%j.log

# Import cost of the proposed default layout against today's layout, with the
# current pbz (layout flags only, no code change). Cold input cache for every
# run (timed_trials.py evicts the input files), output on local disk.
# Submit: sbatch [--nodelist=NODE] workflows/defaults_sweep/import_cost.sh
# Without Slurm: RUN=... WORK=... THREADS=224 bash import_cost.sh

set -euo pipefail
B=/scratch/ucgd/lustre-labs/vollger/users/cade/pbzarr-benchmarking
S=/scratch/ucgd/lustre-labs/vollger/users/cade/defaults-sweep
PY=$B/.pixi/envs/default/bin/python
PBZ=$S/bin/pbz
TIMED=$B/workflows/lib/timed_trials.py
CRAMS=/scratch/ucgd/lustre-labs/vollger/data-sets/1kg-high-coverage/crams
FASTA=/scratch/ucgd/lustre-labs/vollger/assemblies/GRCh38_full_analysis_set_plus_decoy_hla.fa
THREADS=${THREADS:-224}
RUN=${RUN:-$S/import-cost-${SLURM_JOB_ID:-local}}
W=${WORK:-$TMPDIR/import-cost}
OUT=$RUN/results
if [ -e "$RUN" ]; then echo "$RUN exists" >&2; exit 1; fi
mkdir -p "$OUT" "$W"
log() { echo "[$(date +%T)] $*" | tee -a "$OUT/log.txt"; }

declare -A LAYOUT=(
  [old]="--chunk-size 1000000 --column-chunk-size 16"
  [middle]="--chunk-size 65536 --column-chunk-size 16 --shard-size 67108864 --shard-column-size 16"
  [new]="--chunk-size 65536 --column-chunk-size 1 --shard-size 67108864 --shard-column-size 16"
)

# source lists
ls "$B"/results/1kg-sr/bigwig/t64/n32/r1/data/*.bw > "$W/bw32.txt"
ls "$B"/results/1kg-sr/bigwig/t64/n128/r1/data/*.bw > "$W/bw128.txt"
awk -F'\t' -v d="$CRAMS" 'NR > 1 && NR <= 33 {print d "/" $1 ".final.cram"}' "$B/config/1kg-sr.tsv" > "$W/cram32.txt"
{
  echo "start: $(date -Is) on $(hostname), job ${SLURM_JOB_ID:-none}, threads $THREADS"
  echo "pbz: $($PBZ --version) $(git -C $S/pbzarr rev-parse --short HEAD)"
  lscpu | grep -E "Model name|^CPU\(s\)"
  for f in bw32 bw128 cram32; do echo "$f: $(wc -l < "$W/$f.txt") sources, $(xargs du -cb < "$W/$f.txt" | awk 'END {print $1}') bytes"; done
} > "$OUT/context.txt"

run() {
  local input=$1 layout=$2 trials=$3
  local st=$W/$input.$layout.pbz
  local kind=bigwig extra=()
  [[ $input == cram* ]] && { kind=bam; extra=(--fasta "$FASTA"); }
  log "$input $layout ($trials trial(s))"
  mapfile -t data < "$W/$input.txt"
  local dargs=(); for f in "${data[@]}"; do dargs+=(--data "$f"); done
  # shellcheck disable=SC2086
  $PY "$TIMED" "${dargs[@]}" --cache cold --trials "$trials" --before "rm -rf $st" \
    --output "$OUT/$input.$layout.tsv" -- \
    "$PBZ" import $kind "$st" --file-list "$W/$input.txt" --track depth --dim sample \
    -t "$THREADS" --no-progress ${LAYOUT[$layout]} "${extra[@]}" >> "$OUT/runs.log" 2>&1 \
    || log "FAILED $input $layout"
  printf "%s\t%s\t%s\t%s\n" "$input" "$layout" "$(du -sb "$st" 2>/dev/null | cut -f1)" \
    "$(find "$st" -type f 2>/dev/null | wc -l)" >> "$OUT/stores.tsv"
  rm -rf "$st"
}

printf "input\tlayout\tbytes\tfiles\n" > "$OUT/stores.tsv"
for layout in new old middle; do run bw32 "$layout" 2; done
for layout in new old middle; do run cram32 "$layout" 1; done
for layout in new old; do run bw128 "$layout" 1; done

$PY - "$OUT" <<'EOF'
import csv, glob, os, sys
out = sys.argv[1]
rows = []
for f in sorted(glob.glob(f"{out}/*.*.tsv")):
    name = os.path.basename(f)[:-4]
    if name.count(".") != 1:
        continue
    inp, layout = name.split(".")
    for r in csv.DictReader(open(f), delimiter="\t"):
        rows.append({"input": inp, "layout": layout, **r})
with open(f"{out}/imports.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)
EOF
log "done: $OUT"
