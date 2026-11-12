#!/bin/bash
#SBATCH --job-name=pbz-defaults
#SBATCH --partition=reserved
#SBATCH --account=vollger-reserved
# Give the node and size on the command line, e.g.
#   sbatch --nodelist=rw244 --exclusive --mem=1450000M sweep.sh
#   sbatch --nodelist=rw243 -n 1 -c 224 --mem=880G sweep.sh
#SBATCH --time=02:30:00
#SBATCH --output=/scratch/ucgd/lustre-labs/vollger/users/cade/defaults-sweep/logs/sweep-%j.log

# One-hour defaults sweep: re-encode chr20 of two 0.4 stores with zarrs_reencode
# into 0.5 stores of many layouts, then measure size, decode speed, region reads
# and file operations. Submit: sbatch workflows/defaults_sweep/sweep.sh
# Smoke test on a login node: MODE=smoke ROWS=2000000 WORK=/tmp/cade-smoke bash sweep.sh

set -euo pipefail
B=/scratch/ucgd/lustre-labs/vollger/users/cade/pbzarr-benchmarking
S=/scratch/ucgd/lustre-labs/vollger/users/cade/defaults-sweep
D=$B/workflows/defaults_sweep
PY=$B/.pixi/envs/default/bin/python
PY8=$S/pbzarr/.pixi/envs/default/bin/python
PBZ=$S/bin/pbz
RE=$HOME/.cargo/bin/zarrs_reencode
RB=$HOME/.cargo/bin/zarrs_benchmark_read_sync
TIMED=$B/workflows/create/timed_trials.py
MODE=${MODE:-full}
ROWS=${ROWS:-}
RUN=${RUN:-$S/run-${SLURM_JOB_ID:-$MODE}}
W=${WORK:-$TMPDIR/sweep}
L=$RUN/lustre
OUT=$RUN/results
# pin to the CPUs this job may use (a shared node gives a subset)
mapfile -t CPUS < <($PY -c "import os; print(*sorted(os.sched_getaffinity(0)), sep='\n')")
NCPU=${#CPUS[@]}
slot_cpus() { local k=$1; local IFS=,; echo "${CPUS[*]:$((k * 4)):4}"; }
slots=$((NCPU / 4)); [ $slots -gt 64 ] && slots=64
ENC_JOBS=$((NCPU / 8)); [ $ENC_JOBS -gt 48 ] && ENC_JOBS=48; [ $ENC_JOBS -lt 1 ] && ENC_JOBS=1
TRIALS=3; LOCAL_TRIALS=2; THREADS=32
[ "$MODE" = smoke ] && { TRIALS=1; LOCAL_TRIALS=1; ENC_JOBS=4; slots=8; }
declare -A SRC=(
  [dense]=$B/results/1kg-sr/pbz/t16/n128/r1/data/cohort.pbz/depth
  [sparse]=$B/results/1kg-rna-ont/pbz/t64/n60/r1/data/cohort.pbz/depth
)
INPUTS="dense sparse"

if [ -e "$RUN" ]; then echo "$RUN exists; delete it to run again" >&2; exit 1; fi
mkdir -p "$OUT"/{enc,dec,reads} "$L" "$W"
log() { echo "[$(date +%T)] $*" | tee -a "$OUT/log.txt"; }
{
  echo "start: $(date -Is) on $(hostname), job ${SLURM_JOB_ID:-none}, mode $MODE rows ${ROWS:-all}"
  echo "pbz: $($PBZ --version) $(git -C $S/pbzarr rev-parse --short HEAD)"
  echo "zarrs_reencode: $($RE --version)"
  $PY -c "import zarr, xarray, numcodecs; print('bench env: zarr', zarr.__version__, 'xarray', xarray.__version__, 'numcodecs', numcodecs.__version__)"
  $PY8 -c "import pbzarr, zarr; print('pbzarr env:', pbzarr.__version__, 'zarr', zarr.__version__)"
  lscpu | grep -E "Model name|^CPU\(s\)"
  echo "usable cpus: $NCPU, slots $slots, encode jobs $ENC_JOBS"
} > "$OUT/context.txt"

log "extract chr20"
pids=()
for i in $INPUTS; do
  $PY $D/extract_contig.py "${SRC[$i]}" chr20 "$W/in/$i" 32 $ROWS > "$OUT/extract.$i.log" 2>&1 & pids+=($!)
  # the first 2,097,152 positions, for the pinned encode-speed runs
  $PY $D/extract_contig.py "${SRC[$i]}" chr20 "$W/small/$i" 8 2097152 > "$OUT/extract-small.$i.log" 2>&1 & pids+=($!)
done
wait "${pids[@]}"
for i in $INPUTS; do
  len=$($PY -c "import json;print(json.load(open('$W/in/$i/meta.json'))['length'])")
  $PY $D/regions.py chr20 "$len" 2000 7 > "$W/in/$i/regions2000.bed"
  awk 'NR % 10 == 1' "$W/in/$i/regions2000.bed" > "$W/in/$i/regions200.bed"
  n=$($PY -c "import json;print(len(json.load(open('$W/in/$i/meta.json'))['samples']))")
  $PY $D/layouts.py "$n" > "$W/in/$i/layouts.tsv"
  if [ "$MODE" = smoke ]; then
    grep -E "^(codec-bzstd3|codec-pcodec8|ref-control|ref-bench|geo-pcodec8-pc16384-cps4096)	" "$W/in/$i/layouts.tsv" > "$W/in/$i/l.tmp"
    mv "$W/in/$i/l.tmp" "$W/in/$i/layouts.tsv"
  fi
  cp "$W/in/$i/layouts.tsv" "$OUT/layouts.$i.tsv"
done

log "file tree test on Lustre (background)"
(
  $PY $D/filetree.py x 0 0 --header > "$OUT/filetree.tsv"
  if [ "$MODE" = smoke ]; then sizes="1000:1 10000:1"; else sizes="1000:3 10000:3 100000:2 400000:1"; fi
  for nt in $sizes; do $PY $D/filetree.py "$L/filetree" "${nt%:*}" "${nt#*:}" >> "$OUT/filetree.tsv"; done
) > "$OUT/filetree.log" 2>&1 &
FT=$!

# Sets LAYOUT_ARGS (zarrs_reencode layout and codec flags) and GRID for a layout.
layout_args() {
  local row; row=$(awk -F'\t' -v n="$2" '$1 == n' "$W/in/$1/layouts.tsv")
  local chunk shard a2a a2b b2b
  IFS=$'\t' read -r _ GRID chunk shard a2a a2b b2b <<< "$row"
  LAYOUT_ARGS=(-c "$chunk")
  [ "$shard" != - ] && LAYOUT_ARGS+=(-s "$shard")
  [ "$a2a" != - ] && LAYOUT_ARGS+=(--array-to-array-codecs "$a2a")
  [ "$a2b" != - ] && LAYOUT_ARGS+=(--array-to-bytes-codec "$a2b")
  [ "$b2b" != - ] && LAYOUT_ARGS+=(--bytes-to-bytes-codecs "$b2b")
  return 0
}

encode_one() {
  local input=$1 name=$2
  # many jobs run at once; keep each helper's thread pool small
  export RAYON_NUM_THREADS=4
  layout_args "$input" "$name"
  local grid=$GRID
  local st=$W/stores/$input/$name.pbz
  mkdir -p "$W/stores/$input"
  local args=("${LAYOUT_ARGS[@]}" --dimension-names position,sample --concurrent-chunks 8 --validate)
  local log=$OUT/enc/$input.$name.log
  $PY8 $D/make_store.py pre "$W/in/$input/meta.json" "$st" > "$log" 2>&1
  local status=ok
  RAYON_NUM_THREADS=8 /usr/bin/time -f "%e\t%U\t%S\t%M" -o "$OUT/enc/$input.$name.time" \
    $RE "${args[@]}" "$W/in/$input/raw" "$st/depth" >> "$log" 2>&1 || status=reencode_failed
  if [ $status = ok ]; then
    $PY8 $D/make_store.py post "$st" >> "$log" 2>&1 || status=consolidate_failed
    $PBZ validate "$st" >> "$log" 2>&1 || status=invalid
  fi
  local m; m=$($PY $D/measure_store.py "$st" 2>>"$log" || echo "")
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$input" "$name" "$grid" "$status" \
    "$(tail -1 "$OUT/enc/$input.$name.time" 2>/dev/null || printf 'NA\tNA\tNA\tNA')" "$m" > "$OUT/enc/$input.$name.tsv"
}
export -f encode_one layout_args
export W OUT PY PY8 PBZ RE D

log "encode"
for i in $INPUTS; do cut -f1 "$W/in/$i/layouts.tsv" | sed "s/^/$i /"; done \
  | xargs -P "$ENC_JOBS" -n 2 bash -c 'encode_one "$0" "$1"'
{
  printf "input\tname\tgrid\tstatus\tenc_wall_s\tenc_user_s\tenc_sys_s\tenc_max_rss_kb\t"
  $PY $D/measure_store.py x --header
  cat "$OUT"/enc/*.tsv
} > "$OUT/encode.tsv"
log "encode done: $(grep -c "	ok	" "$OUT/encode.tsv" || true) ok of $(($(wc -l < "$OUT/encode.tsv") - 1))"

log "decode speed (4 pinned cores per store, 3 runs)"
ok_store() { awk -F'\t' -v i="$1" -v n="$(basename "$2" .pbz)" '$1 == i && $2 == n && $4 == "ok" {f=1} END {exit !f}' "$OUT/encode.tsv"; }
mapfile -t STORES < <(for i in $INPUTS; do for st in "$W/stores/$i"/*.pbz; do ok_store "$i" "$st" && echo "$st"; done; done)
for ((w = 0; w < ${#STORES[@]}; w += slots)); do
  pids=()
  for ((k = 0; k < slots && w + k < ${#STORES[@]}; k++)); do
    st=${STORES[w + k]}
    (
      for run in 1 2 3; do
        # upstream zarrs fails partial decodes through transpose, so read whole shards there
        inner=--inner-chunks; [[ $st == *ref-bench* ]] && inner=
        line=$(RAYON_NUM_THREADS=4 taskset -c "$(slot_cpus $k)" $RB $inner --concurrent-chunks 4 "$st/depth" 2>&1 | tail -1)
        printf "%s\t%s\t%s\n" "$st" "$run" "$line"
      done > "$OUT/dec/$(basename "$(dirname "$st")").$(basename "$st").tsv"
    ) &
    pids+=($!)
  done
  wait "${pids[@]}"
done
cat "$OUT"/dec/*.tsv > "$OUT/decode.tsv"

log "encode speed (4 pinned cores per layout, 2,097,152 positions, 3 runs)"
mkdir -p "$OUT/spd"
JOBS=()
for i in $INPUTS; do for name in $(cut -f1 "$W/in/$i/layouts.tsv"); do JOBS+=("$i $name"); done; done
for ((w = 0; w < ${#JOBS[@]}; w += slots)); do
  pids=()
  for ((k = 0; k < slots && w + k < ${#JOBS[@]}; k++)); do
    read -r input name <<< "${JOBS[w + k]}"
    (
      layout_args "$input" "$name"
      dst=$W/speed/$input/$name
      for run in 1 2 3; do
        rm -rf "$dst"; mkdir -p "$(dirname "$dst")"
        RAYON_NUM_THREADS=4 /usr/bin/time -f "%e" -o "$dst.time" taskset -c "$(slot_cpus $k)" \
          $RE "${LAYOUT_ARGS[@]}" --concurrent-chunks 4 "$W/small/$input/raw" "$dst" > /dev/null 2>&1 \
          && printf "%s\t%s\t%s\t%s\n" "$input" "$name" "$run" "$(tail -1 "$dst.time")"
      done > "$OUT/spd/$input.$name.tsv"
      rm -rf "$dst"
    ) &
    pids+=($!)
  done
  wait "${pids[@]}"
done
{
  printf "input\tname\trun\twall_s\n"
  cat "$OUT"/spd/*.tsv
} > "$OUT/encode_speed.tsv"
echo "small_cells: $(for i in $INPUTS; do $PY -c "import json;m=json.load(open('$W/small/$i/meta.json'));print('$i', m['length'] * len(m['samples']))"; done | tr '\n' ' ')" >> "$OUT/context.txt"

read_tests() {
  local loc=$1 input=$2 st=$3 trials=$4 tests=$5
  local name; name=$(basename "$st" .pbz)
  local s0; s0=$($PY -c "import json;print(json.load(open('$W/in/$input/meta.json'))['samples'][0])")
  local r2=$W/in/$input/regions2000.bed r1=$W/in/$input/regions200.bed
  for t in $tests; do
    local cmd
    case $t in
      contig_one) cmd=("$PBZ" stat "$st" --track depth -r chr20 -c "$s0" -t $THREADS -o /dev/null --no-header) ;;
      reg_all) cmd=("$PBZ" stat "$st" --track depth -R "$r2" -t $THREADS -o /dev/null --no-header) ;;
      reg_one) cmd=("$PBZ" stat "$st" --track depth -R "$r2" -c "$s0" -t $THREADS -o /dev/null --no-header) ;;
      xr_all) cmd=("$PY" "$D/xr_regions.py" "$st" "$r1" all) ;;
      xr_one) cmd=("$PY" "$D/xr_regions.py" "$st" "$r1" one) ;;
    esac
    $PY "$TIMED" --data "$st" --cache cold --trials "$trials" \
      --output "$OUT/reads/$loc.$input.$name.$t.tsv" -- "${cmd[@]}" >> "$OUT/reads.log" 2>&1 \
      || echo "FAILED $loc $input $name $t" >> "$OUT/reads.log"
  done
}

log "local reads, codec grid"
for i in $INPUTS; do
  for st in "$W/stores/$i"/codec-*.pbz "$W/stores/$i"/ref-*.pbz; do
    [ -e "$st" ] && ok_store "$i" "$st" && read_tests local "$i" "$st" $LOCAL_TRIALS "reg_all reg_one xr_one"
  done
done

log "waiting for the file tree test"
wait $FT || log "file tree test failed (see filetree.log)"

log "copy the geometry stores to Lustre"
for i in $INPUTS; do
  mkdir -p "$L/stores/$i"
  for st in "$W/stores/$i"/geo-*.pbz "$W/stores/$i"/codec-bzstd3.pbz "$W/stores/$i"/codec-pcodec8.pbz "$W/stores/$i"/ref-*.pbz; do
    [ -e "$st" ] && ok_store "$i" "$st" && cp -r "$st" "$L/stores/$i/"
  done
done

log "Lustre reads, cold, one at a time"
for i in $INPUTS; do
  for st in "$L/stores/$i"/*.pbz; do
    read_tests lustre "$i" "$st" $TRIALS "contig_one reg_all reg_one xr_all xr_one"
  done
done

log "collect"
$PY - "$OUT" <<'EOF'
import csv, glob, os, sys
out = sys.argv[1]
rows = []
for f in sorted(glob.glob(f"{out}/reads/*.tsv")):
    loc, inp, name, test = os.path.basename(f)[:-4].split(".", 3)
    for r in csv.DictReader(open(f), delimiter="\t"):
        rows.append({"location": loc, "input": inp, "name": name, "test": test, **r})
if rows:
    with open(f"{out}/reads.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
EOF
echo "end: $(date -Is)" >> "$OUT/context.txt"
log "done: $OUT"
