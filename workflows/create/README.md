# Cohort creation and scan benchmarks

`Snakefile` benchmarks three separate workflows: `pbz.smk`, `d4.smk`, and
`bigwig.smk`. Edit `SAMPLE_COUNTS` and `THREAD_COUNTS` at its top.
`REPLICATES` controls creation trials. `ANALYSIS_REPLICATES` controls scan trials.

Results are written to the git-ignored repository-root `results/` directory,
under `{dataset}/{tool}/t{threads}/n{n_samples}/r{replicate}/`.

Each dataset has a YAML file in the repository-root `config/` directory. Set
its manifest, reference FASTA, reference index, and chromosome sizes. Paths in
YAML and alignment/index paths in the manifest are relative to the repository
root, or absolute. Manifest columns are `sample_id`, `alignment`, and `index`.
Use unique sample IDs that are safe filenames and conventional indexes beside
alignments. Each trial takes the first N manifest rows.

Run from the repository root, selecting a dataset config:

```sh
pixi run snakemake -s workflows/create/Snakefile --cores 64 \
    --resources scan_slot=1 --configfile config/1kg-sr.yaml
```

Each analysis replicate reads creation replicate 1 for the same format,
thread budget, and sample count. For example, `r3/full_scan.tsv` reads
`r1/data/cohort.d4`. Results, logs, and benchmark measurements use the analysis
replicate's directory, so each scan trial has separate outputs.

To run only the analyses, use:

```sh
pixi run snakemake analyses -s workflows/create/Snakefile --cores 64 \
    --resources scan_slot=1 --configfile config/1kg-sr.yaml
```

This target creates replicate 1 data if it is missing. It does not request
other creation replicates or their disk measurements. The default target
still requests all configured creation trials and analysis trials.

The outer `--cores` is the total pool shared by concurrent trials. Each trial
reserves its requested thread count and passes that count to its inner workflow.
For example, `--cores 64` permits four 16-core trials to overlap. Use at least
the largest thread budget. Each trial has a separate working directory.
`--forceall` recreates every inner output, including on retries.

Every timed scan job requests the `scan_slot` resource. Pass
`--resources scan_slot=1` so only one runs at a time. Without the limit,
overlapping jobs share the page cache: one job's cold eviction removes another
job's warm data, and warm reads refill data that a cold job just evicted.

PBZ imports the cohort directly. D4 runs four-thread creation jobs per sample
and merges the tracks with `d4tools merge`. RNA datasets use native sparse
`d4tools create` directly from BAM/CRAM. Other datasets use `mosdepth --d4`.
BigWig uses four-thread Megadepth jobs with `--bigwig --filter-out 1796`.
Both D4 and BigWig creation receive the reference FASTA explicitly.
Inner Snakemake schedules the sample jobs and removes temporary D4 tracks
after the merge.

Sparse D4 defaults to datasets whose name contains `rna` (case insensitive).
It uses `--dict-file d4-zero.dict` for a zero-only dictionary with no primary
table bits, plus `-z` for deflate compression of the secondary table.
Set these dataset config keys to override the defaults:

```yaml
d4_sparse: true
d4_deflate: true  # false disables sparse secondary-table compression
```

Sparse creation uses D4's native alignment coverage and read filters.
Its coverage can differ from mosdepth for spliced reads. `d4_deflate` applies
only to sparse creation.

PBZ and D4 full scans use one BED interval per chromosome from the configured
chromosome sizes. Both formats calculate interval means with their `stat`
commands. D4 disables index shortcuts with `--no-index` for full scans.
All region scans use `regions_bed` from the dataset config.

BigWig scans use `bigwig_scan.py`, the exact pyBigWig helper shared with the
single-sample workflow. Full scans query each chromosome from the BigWig header.
Region scans query each BED interval and preserve their order and duplicates.
Each query uses `stats(type="sum", exact=True)` and divides by the full
interval length. Uncovered bases contribute zero. Exact queries bypass zoom
summaries without expanding every base into a Python array.

Each process scans one sample. Up to the allocated thread count of sample
scans can run concurrently. The helper limits numerical libraries to one
thread before importing pyBigWig, which also loads NumPy. This overrides the
thread limits that Snakemake sets for the entire rule, preventing each worker
from starting its own multithreaded BLAS pool.
Each sample produces `{analysis}/{sample}.bw.tsv`
with headerless `chrom`, `start`, `end`, and `mean` columns, plus a stderr log.
Scan timing includes Python startup, file access, and output writes.
The helper is a workflow input, so edits to it can trigger new scans.

BigWig scan timings now use a different method from the original Megadepth
measurements. Use a fresh dataset result name to keep those measurements
separate. The benchmark paths and aggregate table schema stay the same.

Outer `benchmark.tsv` files include inner Snakemake startup, scheduling, wall
time, and sampled memory. A separate `du -sk` rule records final output size
in KiB in `disk.tsv`. The workflow retains final data, scan outputs, logs,
timings, and disk measurements. It does not measure peak disk use.

D4 retains `benchmarks/encode_sparse_d4/{sample}.tsv` for sparse creation,
or `benchmarks/mosdepth/{sample}.tsv` for dense creation. Both retain
`benchmarks/merge.tsv`.
BigWig retains `benchmarks/megadepth/{sample}.tsv`. Each scan retains
`benchmarks/full_scan.tsv` or `benchmarks/region_scan.tsv`.
The outer `benchmark.tsv` measures total creation time.

Dependencies come from the project Pixi environment. Megadepth is a Linux-only
Pixi dependency, so run the complete cohort workflow on the Linux cluster.
The environment also retains the dependencies for the single-sample workflow.

Aggregate existing benchmark and disk-size TSVs from the repository root:

```sh
pixi run python workflows/create/aggregate_benchmarks.py
```

This writes `results/benchmarks.tsv` and `results/sizes.tsv`, with dataset,
tool, thread budget, sample count, and replicate columns. Benchmark rows also
identify the rule, sample (when applicable), and measurement number. The tables
retain all measurements without averages. Size columns are KiB and GiB.
Analysis result tables are excluded. Pass a different results root
as the optional positional argument.

Plot the aggregate tables with Seaborn, pandas, and Matplotlib:

```sh
pixi run python workflows/create/plot_benchmarks.py
```

This writes creation and analysis figures as PNG and PDF files under
`results/plots/`, separately for each dataset and thread budget. Creation
panels show wall time, peak RSS, and final disk usage. The dashed D4 disk
curve estimates intermediate-plus-final storage as twice final size.
Analysis panels show full-scan and region-scan wall times on a log scale
and peak RSS. Each replicate is drawn separately, without uncertainty bands.

## Creation benchmark

The creation benchmark for the paper is separate from the trials above. It runs
on the exclusive node rw244 and has two parts, each with its own workflows,
configs and batch scripts. The design is in
`docs/superpowers/specs/2026-10-20-creation-benchmark-design.md`.

- `workflows/single` is part 1: one depth track from one CRAM, for PBZ,
  mosdepth (D4) and Megadepth (BigWig), at several thread counts. Config:
  `config/single.yaml`.
- `workflows/cohort` is part 2: one cohort store built from the first N
  manifest samples, with one workflow per tool (`pbz.smk`, `d4.smk`,
  `bigwig.smk`), plus the cold scans. Config: `config/cohort.yaml`.
- `workflows/lib/evict.py` drops files from the page cache. Both parts use it.

Run everything from the cluster checkout,
`/scratch/ucgd/lustre-labs/vollger/users/cade/pbzarr-benchmarking`. Slurm does
not create the log folders, so create them once:

```sh
mkdir -p results/single/logs results/single-smoke/logs results/cohort/logs results/cohort-smoke/logs
```

Copy the files to the cluster from this repo. Do not use `--delete`: the cluster
checkout has its own `config/1kg-sr.tsv`, `codecs.json`, `random.bed` and
`results/`.

```sh
rsync -av --relative workflows/lib workflows/single workflows/cohort \
    config/single.yaml config/single-smoke.yaml config/cohort.yaml \
    config/cohort-smoke.yaml config/packing-smoke.yaml \
    redwood:/scratch/ucgd/lustre-labs/vollger/users/cade/pbzarr-benchmarking/
```

Check a workflow on a login node with a dry run. Part 1 expects 153 `run_*`
jobs:

```sh
export PATH=$PWD/.pixi/envs/default/bin:$PATH
snakemake -s workflows/single/Snakefile --configfile config/single.yaml \
    --directory /tmp/single-dry --cores 256 --resources bench=1 -n --printshellcmds
```

For part 2, `/tmp/s32.txt` holds the first 32 sample IDs of the manifest and
`/tmp/ref.fa` is any placeholder file. Use `pbz`, `d4` or `bigwig` for `<tool>`:

```sh
snakemake -s workflows/cohort/<tool>.smk build scans --configfile config/cohort.yaml \
    --config n_samples=32 samples_file=/tmp/s32.txt results=/tmp/r build=1 \
    threads_per_job=2 reference_local=/tmp/ref.fa --directory /tmp/<tool>-dry \
    --cores 256 --resources pack=100 bench=1 -n --printshellcmds
```

**Smoke tests.** Each needs rw244 to itself, so make sure no other job runs there
first. Part 1 takes about 20 minutes, part 2 more than an hour. Submit again after
a smoke run and nothing is redone: part 1 skips finished runs and part 2 prints
`already done`.

```sh
CONFIG=config/single-smoke.yaml sbatch workflows/single/part1.sbatch
CONFIG=config/cohort-smoke.yaml PACKING=config/packing-smoke.yaml \
    sbatch workflows/cohort/part2.sbatch 32 1
```

**Real runs.**

1. Part 1: `sbatch workflows/single/part1.sbatch`. If the job dies, submit it again.
   It writes `results/single/runs.tsv` (one row per run), `samples.tsv` and
   `context.txt`.
2. Choose the D4 and BigWig packing from part 1's medians:
   `python workflows/cohort/packing.py choose results/single/runs.tsv > config/packing.yaml`.
3. Part 2: `workflows/cohort/submit_part2.sh` submits every cohort build as its
   own job, one after another: N = 32, 64, 128, 256, 512, 1024, then builds 2
   and 3 at N = 256. A job that died starts over when submitted again; a job
   with a `DONE` file exits. Each job writes
   `results/cohort/n<N>/b<build>/` with `builds.tsv`, `merge.tsv` (D4 only),
   `scans/`, `context.txt`, `samples.txt`, `order.txt`, `logs/` and `DONE`.

Part 1's batch job asks Slurm for 1450000M of memory, not `--mem=0`: all of
rw244's memory is above the account's group limit, so `--mem=0` is rejected.
`--exclusive` still gives the job the whole node.
