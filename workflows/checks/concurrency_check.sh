#!/usr/bin/env bash
# Do different tools slow each other down when they run at the same time?
# Throwaway check. Each tool runs alone, then all three run at once, twice, at
# 64 threads. Each tool reads its own CRAM, so no tool warms the page cache for
# another. The CRAMs are cleared from the cache before every round.
# Submit from the repo root:  sbatch workflows/checks/concurrency_check.sh
#SBATCH --job-name=concurrency-check
#SBATCH --partition=reserved
#SBATCH --account=vollger-reserved
#SBATCH --nodelist=rw244
#SBATCH --exclusive
#SBATCH --mem=1450000M
#SBATCH --time=04:00:00
#SBATCH --output=results/checks/concurrency-check-%j.log

set -euo pipefail

R=/scratch/ucgd/lustre-labs/vollger/users/cade/pbzarr-benchmarking
CRAMS=/scratch/ucgd/lustre-labs/vollger/data-sets/1kg-high-coverage/crams
REF=/scratch/ucgd/lustre-labs/vollger/assemblies/GRCh38_full_analysis_set_plus_decoy_hla.fa
PBZ=/uufs/chpc.utah.edu/common/HIPAA/u1364833/pbzarr/target/release/pbz
export PATH=$R/.pixi/envs/default/bin:$PATH
THREADS=64
declare -A SAMPLE=([pbz]=HG00475 [mosdepth]=HG01748 [megadepth]=HG00244)

OUT=$R/results/checks/concurrency/$SLURM_JOB_ID
W=$TMPDIR/concurrency
mkdir -p "$OUT" "$W/out"
cp "$REF" "$REF.fai" "$W/"
ref=$W/$(basename "$REF")
printf 'mode\trep\ttool\tsample\tthreads\twall_s\tuser_s\tsys_s\tmax_rss_kb\n' > "$OUT/runs.tsv"

command_for() {  # command_for TOOL: the tool's command, writing under $W/out/TOOL
    local tool=$1 sample=${SAMPLE[$1]} out=$W/out/$1
    local cram=$CRAMS/$sample.final.cram
    case $tool in
        pbz) echo "$PBZ import bam $cram:$sample --reference $ref --threads $THREADS --track depth --codecs @$R/n1codecs.json --no-progress -o $out/s.pbz" ;;
        mosdepth) echo "mosdepth -t $THREADS --fasta $ref --d4 $out/s $cram" ;;
        megadepth) echo "megadepth $cram --fasta $ref --bigwig --filter-out 1796 --threads $THREADS --prefix $out/s" ;;
    esac
}

start() {  # start TOOL MODE REP: one timed run in the background
    local tool=$1 mode=$2 rep=$3
    rm -rf "$W/out/$tool"
    mkdir -p "$W/out/$tool"
    (
        /usr/bin/time -f '%e\t%U\t%S\t%M' -o "$W/out/$tool.time" \
            sh -c "$(command_for "$tool") && sync" > "$W/out/$tool.log" 2>&1
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$mode" "$rep" "$tool" "${SAMPLE[$tool]}" "$THREADS" \
            "$(cat "$W/out/$tool.time")" >> "$OUT/runs.tsv"
    ) &
}

clear_cache() {
    local paths=()
    for tool in pbz mosdepth megadepth; do
        paths+=("$CRAMS/${SAMPLE[$tool]}.final.cram" "$CRAMS/${SAMPLE[$tool]}.final.cram.crai")
    done
    python "$R/workflows/lib/evict.py" "${paths[@]}"
}

for rep in 1 2; do
    for tool in pbz mosdepth megadepth; do
        clear_cache
        start "$tool" alone "$rep"
        wait
    done
    clear_cache
    for tool in pbz mosdepth megadepth; do
        start "$tool" together "$rep"
    done
    wait
done

column -t "$OUT/runs.tsv"
