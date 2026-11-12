#!/bin/bash
# Build today's layout ([1,000,000 x 16], blosc zstd5 shuffle, no shard) with
# zarr-python, because zarrs_reencode keeps the input chunks when -s is not
# given. Then run the same measurements as sweep.sh on it, as "ref-today".
# usage: RUN=... WORK=... bash control.sh   (after sweep.sh, same RUN and WORK)

set -euo pipefail
B=/scratch/ucgd/lustre-labs/vollger/users/cade/pbzarr-benchmarking
S=/scratch/ucgd/lustre-labs/vollger/users/cade/defaults-sweep
D=$B/workflows/defaults_sweep
PY=$B/.pixi/envs/default/bin/python
PY8=$S/pbzarr/.pixi/envs/default/bin/python
PBZ=$S/bin/pbz
RB=$HOME/.cargo/bin/zarrs_benchmark_read_sync
TIMED=$B/workflows/create/timed_trials.py
RUN=${RUN:?}; W=${WORK:?}; OUT=$RUN/results; L=$RUN/lustre
TRIALS=3; LOCAL_TRIALS=2; THREADS=32
NAME=ref-today
export RAYON_NUM_THREADS=8
log() { echo "[$(date +%T)] $*" | tee -a "$OUT/log.txt"; }

for i in dense sparse; do
  st=$W/stores/$i/$NAME.pbz
  rm -rf "$st"
  log "control $i: write"
  $PY8 $D/make_store.py pre "$W/in/$i/meta.json" "$st"
  t0=$(date +%s.%N)
  $PY - "$W/in/$i/raw" "$st" <<'EOF'
import sys
from multiprocessing import Pool
import zarr
from zarr.codecs import BloscCodec
raw_path, st = sys.argv[1], sys.argv[2]
raw = zarr.open_array(raw_path, mode="r")
n_pos, n = raw.shape
zarr.create_array(st, name="depth", shape=raw.shape, chunks=(1_000_000, min(16, n)), dtype="int32",
                  fill_value=0, dimension_names=["position", "sample"],
                  compressors=BloscCodec(cname="zstd", clevel=5, shuffle="shuffle", typesize=4))
def copy(a):
    b = min(a + 1_000_000, n_pos)
    zarr.open_array(f"{st}/depth", mode="r+")[a:b] = zarr.open_array(raw_path, mode="r")[a:b]
with Pool(32) as p:
    list(p.imap_unordered(copy, range(0, n_pos, 1_000_000)))
EOF
  wall=$(echo "$(date +%s.%N) - $t0" | bc)
  $PY8 $D/make_store.py post "$st"
  status=ok; $PBZ validate "$st" >> "$OUT/control.$i.log" 2>&1 || status=invalid
  printf "%s\t%s\tcodec\t%s\t%s\tNA\tNA\tNA\t%s\n" "$i" "$NAME" "$status" "$wall" "$($PY $D/measure_store.py "$st")" >> "$OUT/encode.tsv"

  log "control $i: decode"
  for run in 1 2 3; do
    line=$(RAYON_NUM_THREADS=4 taskset -c 0-3 $RB --concurrent-chunks 4 "$st/depth" 2>&1 | tail -1)
    printf "%s\t%s\t%s\n" "$st" "$run" "$line"
  done >> "$OUT/decode.tsv"

  s0=$($PY -c "import json;print(json.load(open('$W/in/$i/meta.json'))['samples'][0])")
  r2=$W/in/$i/regions2000.bed; r1=$W/in/$i/regions200.bed
  mkdir -p "$L/stores/$i"; rm -rf "$L/stores/$i/$NAME.pbz"; cp -r "$st" "$L/stores/$i/"
  for loc in local lustre; do
    if [ $loc = local ]; then ls=$st; trials=$LOCAL_TRIALS; tests="reg_all reg_one xr_one"
    else ls=$L/stores/$i/$NAME.pbz; trials=$TRIALS; tests="contig_one reg_all reg_one xr_all xr_one"; fi
    log "control $i: $loc reads"
    for t in $tests; do
      case $t in
        contig_one) cmd=("$PBZ" stat "$ls" --track depth -r chr20 -c "$s0" -t $THREADS -o /dev/null --no-header) ;;
        reg_all) cmd=("$PBZ" stat "$ls" --track depth -R "$r2" -t $THREADS -o /dev/null --no-header) ;;
        reg_one) cmd=("$PBZ" stat "$ls" --track depth -R "$r2" -c "$s0" -t $THREADS -o /dev/null --no-header) ;;
        xr_all) cmd=("$PY" "$D/xr_regions.py" "$ls" "$r1" all) ;;
        xr_one) cmd=("$PY" "$D/xr_regions.py" "$ls" "$r1" one) ;;
      esac
      $PY "$TIMED" --data "$ls" --cache cold --trials "$trials" \
        --output "$OUT/reads/$loc.$i.$NAME.$t.tsv" -- "${cmd[@]}" >> "$OUT/reads.log" 2>&1 \
        || echo "FAILED $loc $i $NAME $t" >> "$OUT/reads.log"
    done
  done
done

log "control: collect"
$PY - "$OUT" <<'EOF'
import csv, glob, os, sys
out = sys.argv[1]
rows = []
for f in sorted(glob.glob(f"{out}/reads/*.tsv")):
    loc, inp, name, test = os.path.basename(f)[:-4].split(".", 3)
    for r in csv.DictReader(open(f), delimiter="\t"):
        rows.append({"location": loc, "input": inp, "name": name, "test": test, **r})
with open(f"{out}/reads.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)
EOF
log "control done"
