#!/usr/bin/env bash
# Where do the tools' per-base depths differ? Throwaway check.
# Runs each tool (and a few variants of mosdepth and d4tools) on the chr22 reads
# of one alignment file, then compares chr22 base by base with depth_compare.py.
# Submit from the repo root:
#   sbatch workflows/checks/depth_compare.sh                            (the smoke CRAM)
#   sbatch workflows/checks/depth_compare.sh SAMPLE ALIGNMENT REFERENCE
# Give REFERENCE as '' for a BAM: the tools then get no reference option.
#SBATCH --job-name=depth-compare
#SBATCH --partition=reserved
#SBATCH --account=vollger-reserved
#SBATCH --nodelist=rw244
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=results/checks/depth-compare-%j.log

set -euo pipefail

R=/scratch/ucgd/lustre-labs/vollger/users/cade/pbzarr-benchmarking
SAMPLE=${1:-HG00096}
ALIGNMENT=${2:-$R/smoke-crams/$SAMPLE.final.cram}
REF=${3-/scratch/ucgd/lustre-labs/vollger/assemblies/GRCh38_full_analysis_set_plus_decoy_hla.fa}
PBZ=/uufs/chpc.utah.edu/common/HIPAA/u1364833/pbzarr/target/release/pbz
# samtools is not in the pixi env; this is CHPC's samtools/1.16 module.
SAMTOOLS=/uufs/chpc.utah.edu/sys/installdir/samtools/1.16/bin/samtools
export PATH=$R/.pixi/envs/default/bin:$PATH
T=4

OUT=$R/results/checks/depth-compare/$SLURM_JOB_ID
W=$TMPDIR/depth-compare
mkdir -p "$OUT" "$W"

# The reference options, empty for a BAM.
fasta=() d4_ref=() pbz_ref=()
if [ -n "$REF" ]; then
    fasta=(--fasta "$REF") d4_ref=(-r "$REF") pbz_ref=(--reference "$REF")
fi

# The chr22 reads only, so every tool runs on a small file.
if [ -n "$REF" ]; then
    CRAM=$W/$SAMPLE.chr22.cram
    $SAMTOOLS view -@ $T -C -T "$REF" -o "$CRAM" "$ALIGNMENT" chr22
else
    CRAM=$W/$SAMPLE.chr22.bam
    $SAMTOOLS view -@ $T -b -o "$CRAM" "$ALIGNMENT" chr22
fi
$SAMTOOLS index "$CRAM"

run() {  # run LABEL CMD...: run one tool, log its output, record its wall time
    local label=$1; shift
    echo "== $label: $*"
    /usr/bin/time -f "$label\t%e s\t%M KB" -a -o "$OUT/times.tsv" "$@" > "$W/$label.log" 2>&1 \
        || { echo "$label failed"; tail -20 "$W/$label.log"; }
}

# The pairs that should match exactly, from reading the sources:
#   d4tools (defaults: MAPQ >= 60, no flag filter) = mosdepth -x -Q 60 -F 0
#   d4tools -q 0 -F '~1796'                      = mosdepth -x (deletions and mate overlaps counted)
#   pbz --overlap none --count-deletions         = mosdepth -x
#   pbz (defaults)                               = mosdepth (defaults)
run mosdepth          mosdepth -t $T "${fasta[@]}" --d4 $W/mosdepth $CRAM
run mosdepth_x        mosdepth -t $T "${fasta[@]}" --d4 -x $W/mosdepth_x $CRAM
run mosdepth_x_Q60_F0 mosdepth -t $T "${fasta[@]}" --d4 -x -Q 60 -F 0 $W/mosdepth_x_Q60_F0 $CRAM
run d4tools           d4tools create -t $T "${d4_ref[@]}" $CRAM $W/d4tools.d4
run d4tools_q0_F1796  d4tools create -t $T "${d4_ref[@]}" -q 0 -F '~1796' $CRAM $W/d4tools_q0_F1796.d4
run megadepth         megadepth $CRAM "${fasta[@]}" --bigwig --filter-out 1796 --threads $T --prefix $W/megadepth
run pbz               $PBZ import bam $CRAM:$SAMPLE "${pbz_ref[@]}" --threads $T --track depth \
                          --codecs @$R/n1codecs.json --no-progress -o $W/pbz.pbz
run pbz_nooverlap_del $PBZ import bam $CRAM:$SAMPLE "${pbz_ref[@]}" --threads $T --track depth \
                          --overlap none --count-deletions \
                          --codecs @$R/n1codecs.json --no-progress -o $W/pbz_nooverlap_del.pbz

# depth_compare.py reads <label>.d4 / .bw / .pbz
for f in $W/*.per-base.d4; do mv "$f" "${f%.per-base.d4}.d4"; done
mv $W/megadepth.all.bw $W/megadepth.bw

python $R/workflows/checks/depth_compare.py $W chr22 mosdepth | tee "$OUT/compare.txt"
cp $W/*.bed $W/*.log "$OUT/"
