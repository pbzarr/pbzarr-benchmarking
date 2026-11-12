# Cohort workflow for D4: one mosdepth run per sample, then one `d4tools merge`,
# then the scans.
#   snakemake -s workflows/cohort/d4.smk build   (timed by the driver, with timed_trials.py)
#   snakemake -s workflows/cohort/d4.smk scans   (each scan timed by timed_trials.py)
# See common.smk.

include: "common.smk"

wildcard_constraints:
    sample="[^/]+"

EVICT = str(ROOT / "workflows/lib/evict.py")
# The cohort copy of workflows/create/over10_d4.py: window and dictionary sized to N.
OVER10_D4 = str(ROOT / "workflows/cohort/over10_d4.py")
PER_SAMPLE = [f"tmp/{sample}.per-base.d4" for sample in SAMPLES]


rule build:
    input:
        D4_FILE


rule scans:
    input:
        expand(RESULTS + "/scans/d4.{analysis}.tsv", analysis=["full_scan", "region_scan", "over10"])


# Each job declares threads_per_job threads. Each sample is declared as an
# output, so a missing file fails the build.
rule mosdepth:
    input:
        alignment=lambda wildcards: alignment(wildcards.sample),
        index=lambda wildcards: index(wildcards.sample)
    output:
        temp("tmp/{sample}.per-base.d4")
    threads: int(config["threads_per_job"])
    params:
        reference=reference_option("--fasta"),
        prefix="tmp/{sample}"
    shell:
        "mosdepth -t {threads} {params.reference} --d4 {params.prefix} {input.alignment:q}"


# D4 is the only tool that writes files during the build and reads them back.
# Before the merge, the per-sample files are flushed to disk and cleared from the
# page cache, so the merge reads them cold at every N: at small N they would
# still be cached, and at N = 1024 (2.7 TB) they cannot be.
# The merge is timed on its own with /usr/bin/time (wall s, user s, sys s, max
# RSS KB), with `sync` inside the timed command. Paths have no spaces, so they
# are not quoted inside sh -c '...'. merge.tsv also gets the per-sample bytes:
# they sit on disk next to the merged file, so D4's peak disk is the two added.
rule merge:
    input:
        PER_SAMPLE
    output:
        D4_FILE
    params:
        sources=[f"tmp/{sample}.per-base.d4:{sample}" for sample in SAMPLES],
        evict=EVICT,
        results=RESULTS,
        n=config["n_samples"]
    shell:
        r"""
        per_sample_bytes=$(du -cb {input:q} | tail -n 1 | cut -f1)
        sync
        python {params.evict:q} {input:q}
        /usr/bin/time -f '%e\t%U\t%S\t%M' -o merge.time sh -c 'd4tools merge {params.sources} {output} && sync'
        read wall user sys rss < merge.time
        bytes=$(du -sb {output:q} | cut -f1)
        printf 'n\twall_s\tuser_s\tsys_s\tmax_rss_kb\tout_bytes\tper_sample_bytes\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            {params.n} "$wall" "$user" "$sys" "$rss" "$bytes" "$per_sample_bytes" > {params.results:q}/merge.tsv
        """


# The scans below run the rules d4_stat and d4_over10 of
# workflows/create/analysis.smk through the same timed_trials.py. The over10
# script is the cohort copy (workflows/cohort/over10_d4.py, sized to N). What
# differs: the data is in the working directory, threads are CORES, the cache is
# always cold, and trials come from the config.
#
# Everything a scan writes (stat tables, over10 outputs, logs) goes in the working
# directory on local disk. Only the timing table goes to
# {results}/scans/d4.{analysis}.tsv, written by timed_trials.py after the last
# trial, so no Lustre write lands inside a timed scan. Scans run one at a time
# (resources: bench=1), because overlapping jobs evict and refill each other's
# page cache.

rule d4_stat:
    input:
        d4=D4_FILE,
        regions=lambda wc: "full_scan.bed" if wc.analysis == "full_scan" else REGIONS_BED,
        timer=TIMED_TRIALS
    output:
        tsv="scan/d4.{analysis}.tsv",
        bench=RESULTS + "/scans/d4.{analysis}.tsv"
    threads: CORES
    resources:
        bench=1
    params:
        python=PYTHON,
        trials=TRIALS,
        scan_flag=lambda wc: "--no-index" if wc.analysis == "full_scan" else ""
    log:
        "logs/d4.{analysis}.log"
    shell:
        """
        # d4tools stat prints to stdout, so sh -c sends it to the file. Paths have
        # no spaces, so they are not quoted inside sh -c '...'.
        {params.python:q} {input.timer:q} --data {input.d4:q} --cache cold \
            --trials {params.trials} --output {output.bench:q} -- \
            sh -c 'd4tools stat {params.scan_flag} -s mean -H -t {threads} -r {input.regions} {input.d4} > {output.tsv}' \
            2> {log:q}
        """


# Per-base fraction of samples with depth >= 10x (gnomAD's over_10), written as a
# D4 file by a plain script that uses a multiprocessing pool.
# --before removes the previous trial's output outside the timed window.
rule d4_over10:
    input:
        d4=D4_FILE,
        script=OVER10_D4,
        timer=TIMED_TRIALS
    output:
        out="scan/d4.over10.d4",
        bench=RESULTS + "/scans/d4.over10.tsv"
    threads: CORES
    resources:
        bench=1
    params:
        python=PYTHON,
        trials=TRIALS
    log:
        "logs/d4.over10.log"
    shell:
        """
        # The script runs {threads} worker processes. pyd4 also starts a pool as
        # large as the node (512 threads) in each one, and 256 of them hit the
        # per-user thread limit and hang. One thread each leaves the parallelism
        # to the worker processes.
        RAYON_NUM_THREADS=1 \
        {params.python:q} {input.timer:q} --data {input.d4:q} --cache cold \
            --trials {params.trials} --output {output.bench:q} --before 'rm -rf {output.out}' -- \
            {params.python:q} {input.script:q} {input.d4:q} {output.out:q} {threads} \
            2> {log:q}
        """
