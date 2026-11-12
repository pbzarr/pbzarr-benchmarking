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

- **Part 1** (`workflows/single/`) times making one depth track from one
  alignment file, for PBZ, mosdepth (D4), d4tools create (D4) and Megadepth
  (BigWig), at several thread counts, then checks where the tools' depths
  agree. One config per data type: `config/single.yaml` (short read),
  `config/single-ont.yaml` and `config/single-rna.yaml`. Submit one job per
  config with `sbatch workflows/single/part1.sbatch CONFIG`.
- **Part 2** (`workflows/cohort/`, config `config/cohort.yaml`) builds one
  cohort store from the first N manifest samples, once for each tool, and runs
  the cold scans on it. There is one workflow per tool (`pbz.smk`, `d4.smk`,
  `bigwig.smk`, which share `common.smk`), each with a `build` target and a
  `scans` target. One config per dataset (`config/cohort*.yaml`). Submit one
  job with `sbatch workflows/cohort/part2.sbatch CONFIG N`, or one job for each
  of the config's `sample_counts` with `workflows/cohort/submit_part2.sh CONFIG`.
- `workflows/lib/evict.py` drops files from the page cache. Both parts use it.

Run everything from the cluster checkout,
`/scratch/ucgd/lustre-labs/vollger/users/cade/pbzarr-benchmarking`. Slurm does
not create the log folders, so create them once:

```sh
mkdir -p results/logs
```

### How things are measured

- **Part 1:** each run is one command inside `/usr/bin/time -f '%e\t%U\t%S\t%M'`,
  with `sync` in the timed command. The alignment file and its index are dropped
  from the page cache first. The row in `runs.tsv` has wall, user and system
  seconds, peak RSS in KB, and the size and file count of the output.
- **Part 1 correctness check:** after the timed runs, each tool runs once more,
  untimed, at its largest thread count. `workflows/single/agreement.py`
  compares the outputs base by base on chr22, for every pair of tools.
  `d4tools create` gets the other tools' read filters (`-q 0 -F '~1796'`); for
  RNA it also gets `--sparse`, as in the D4 paper. It still counts deletions,
  and on RNA the spliced-out introns, so it differs from the others there.
- **Part 2 builds:** `workflows/create/timed_trials.py` times the whole
  `snakemake build` call (plus `sync`), one trial, with every CRAM and `.crai`
  of the cohort dropped from the page cache first. Wall time includes
  Snakemake's start-up. Peak RSS is summed over the process tree and sampled.
- **Part 2 scans:** `timed_trials.py` too, with a cold cache before every trial
  (`scan_trials` trials), one scan at a time.
- **D4 merge:** `d4tools merge` inside `/usr/bin/time`, in `d4.smk`. It goes to
  `merge.tsv`. The per-sample files are dropped from the page cache before it:
  D4 is the only tool that writes files in the build and reads them back.
- **Peak disk:** computed from file sizes, not sampled. It is the output size,
  plus for D4 the size of the per-sample files that sit next to the merged
  file.

### Thread counts

PBZ always uses all `cores` (256) in one process. D4 and BigWig run one
mosdepth or Megadepth process per sample, all at once, with `cores / N` threads
each (rounded down). From N = `cores` up, each gets 1 thread and Snakemake runs
`cores` of them at a time. `snakemake` always gets `--cores 256`. Details are in
`workflows/cohort/common.smk`.

### Copy to the cluster

Do not use `--delete`: the cluster checkout has its own `config/1kg-sr.tsv`,
`codecs.json` (cohort, 2D), `n1codecs.json` (part 1, 1D single-sample store),
`random.bed` and `results/`.

Part 2 also reads scripts from `workflows/create` (`timed_trials.py`,
`bigwig_scan.py`, `over10_pbz.py`, `proposed_api/`; `over10_d4.py` is the
original of the copy in `workflows/cohort`), so the command names them one by one.
Never sync the whole `workflows/create` folder: the cluster copy has files that
exist only there (`disk.smk`, an edited `pbz.smk`), and a folder sync would
overwrite them.

```sh
rsync -av --relative workflows/lib workflows/single workflows/cohort \
    workflows/create/timed_trials.py workflows/create/bigwig_scan.py \
    workflows/create/over10_pbz.py workflows/create/over10_d4.py \
    workflows/create/proposed_api \
    config/single*.yaml config/smoke-manifest.tsv \
    config/cohort*.yaml \
    redwood:/scratch/ucgd/lustre-labs/vollger/users/cade/pbzarr-benchmarking/
```

### Dry runs

Check a workflow on a login node. A new part 1 config expects 48 `run_*` jobs
(4 tools, 4 thread counts, 3 repeats; 4 with `config/single-smoke.yaml`), 4
`check_*` jobs and 1 `agreement` job. Finished runs are left out, so a config
that already ran shows only the runs still to do:

```sh
export PATH=$PWD/.pixi/envs/default/bin:$PATH
snakemake -s workflows/single/Snakefile --configfile config/single.yaml \
    --directory /tmp/single-dry --cores 256 --resources bench=1 -n
```

For part 2, `/tmp/samples.txt` holds the sample IDs (as many as `n_samples`, one
per line) and `/tmp/ref.fa` is any placeholder file. Use `pbz`, `d4` or `bigwig`
for `<tool>`:

```sh
snakemake -s workflows/cohort/<tool>.smk build scans --configfile config/cohort-ont.yaml \
    --config n_samples=4 samples_file=/tmp/samples.txt results=/tmp/r \
    threads_per_job=64 reference_local=/tmp/ref.fa --directory /tmp/<tool>-dry \
    --cores 256 --resources bench=1 -n --printshellcmds
```

### Part 1 smoke test

It needs rw244 to itself, so make sure no other job runs there first. It uses
four chr22-only CRAMs in `smoke-crams/` (listed in `config/smoke-manifest.tsv`)
and takes a few minutes. Submit again after a smoke run and only the check is redone.

```sh
sbatch workflows/single/part1.sbatch config/single-smoke.yaml
```

### Real runs

1. Part 1: `sbatch workflows/single/part1.sbatch CONFIG`, once for each of
   `config/single.yaml`, `config/single-ont.yaml` and `config/single-rna.yaml`.
   If a job dies, submit it again: runs that already have a result file are
   skipped.
2. Part 2: `workflows/cohort/submit_part2.sh CONFIG` submits one job for each
   of the config's `sample_counts`. Each job may run 3 days. A job stops at once
   if its results folder exists: delete the folder of a job that died, then
   submit again. Each tool's scans stop after 16 hours; if they fail or time
   out, the job notes it in `failed.txt` and goes on, so the build rows survive.

### Where results land

Part 1, under each config's `results` folder (`results/single/`,
`results/single-ont/`, `results/single-rna/`, `results/single-smoke/`). The job
logs go to `results/logs/part1-<job id>.log`.

- `runs.tsv`: one row per run (tool, sample, threads, repeat, input bytes, wall,
  user and system seconds, peak RSS KB, output bytes and files). The input
  bytes column is called `cram_bytes` for every data type, BAM input included,
  so the tables line up.
- `runs/<tool>/<sample>/t<threads>/r<repeat>.tsv`: the one-row file for each
  run, which `runs.tsv` joins.
- `samples.tsv`: the samples picked, by alignment file size.
- `agreement/<sample>.tsv`: the correctness check, one row per pair of tools:
  chr22 bases covered by either tool, how many of them differ, the percentage,
  and the largest difference. The files of all data types stack into one table.
- `context.txt`: each job adds its start time, host, job id, pbz version and
  commit, and tool versions. The jobs before 2026-10-27 wrote more lines.
- `config.yaml`: a copy of the config file used. If the file changes before a
  later submit, the older copy is kept as `config.yaml.~1~`.

Part 2, one folder for each job, `<results>/n<N>/` (the runs before
2026-10-26 used `n<N>/b<build>/`):

- `builds.tsv`: one row per tool, in build order: threads per job, bytes read,
  wall seconds, peak RSS in MiB, MiB read from storage (`io_in`), CPU seconds,
  and output bytes and files.
- `merge.tsv`: the D4 merge, one row. `per_sample_bytes` is the size of the
  per-sample D4 files; D4's peak disk is that plus the merged file.
- `scans/<tool>.<analysis>.tsv`: the timing tables of the scans.
- `failed.txt`: only if a tool's scans failed or timed out (exit status 124).
  The build rows are still valid; that tool's scan tables may be missing.
- `timing/build-<tool>.tsv`: the raw `timed_trials.py` row of each build.
- `samples.txt`, `context.txt` (start time, host, job id, pbz version and
  commit, tool versions), `config.yaml` (a copy of the config file used),
  `logs/` and `DONE` (the end time).

The analysis must read only folders that have `DONE`.

### Notes for the methods

- **RSS includes file-backed pages.** `max_rss` and `max_rss_mib` count every
  resident page of the process, including pages of files that the tool maps
  into memory. This matters for D4's writer: `mosdepth --d4` maps its output
  file, so its RSS is about the size of that file, whether or not the process
  allocated that memory. Measured on one smoke CRAM (HG00096, 4 threads, a
  `/proc/<pid>/status` sample each second, on a `ucgd-rw` node): `mosdepth --d4`
  peaked at 3.06 GiB VmRSS, of which 0.43 GiB was RssAnon and 2.63 GiB RssFile
  (its output file is 2,818,589,910 bytes, the same size). `megadepth --bigwig`
  peaked at 0.41 GiB VmRSS, with 0.40 GiB RssAnon and 5.8 MiB RssFile. Say so
  wherever D4's peak memory is compared with the other tools.
- **D4 over10 uses a window and a dictionary width sized to N.**
  `workflows/cohort/over10_d4.py` reads windows of `min(1 Mb, 128 Mb / N)`
  bases, so a worker holds about the memory it holds at N = 128, and it sets the
  dictionary bits to `ceil(log2(N + 1))`, so counts 0 to N stay in the primary
  table. Up to N = 128 the window is the 1 Mb of the older runs.
- **D4 region scans use no index**, as in the D4 paper's benchmarks. The build
  never runs `d4tools index`.
- **D4's per-sample phase (build minus merge) includes eviction and temp-file
  deletion.** Besides the mosdepth runs it holds the `sync` and page-cache
  eviction before the merge, and deleting the per-sample files afterwards.
- **Kernel writeback and encryption CPU is not in `cpu_s`.** `cpu_s` is the user
  and system time of the tool's own processes; kernel threads that write dirty
  pages back or encrypt them are not counted.
- **A `mean_load` above 100 x cores means hyperthreads were used.** `mean_load`
  in the timing files is `100 * cpu_time / s`, so one fully busy CPU is 100. The
  node has 512 logical CPUs; with `cores` at 256, a value above 25600 means the
  process ran on the second thread of some cores.

Both batch jobs ask Slurm for 1450000M of memory, not `--mem=0`: all of rw244's
memory is above the account's group limit, so `--mem=0` is rejected.
`--exclusive` still gives the job the whole node.
