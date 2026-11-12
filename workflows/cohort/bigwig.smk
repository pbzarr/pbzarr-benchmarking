# Cohort workflow for BigWig: one Megadepth run per sample, then the scans.
#   snakemake -s workflows/cohort/bigwig.smk build   (timed by the driver, with timed_trials.py)
#   snakemake -s workflows/cohort/bigwig.smk scans   (each scan timed by timed_trials.py)
# See common.smk.

include: "common.smk"

wildcard_constraints:
    sample="[^/]+"

BIGWIG_SCAN = str(CREATE / "bigwig_scan.py")
OVER10_BIGWIG = str(CREATE / "over10_bigwig.py")


rule build:
    input:
        BIGWIGS


rule scans:
    input:
        expand(RESULTS + "/scans/bigwig.{analysis}.tsv", analysis=["full_scan", "region_scan", "over10"])


# Each job declares threads_per_job threads. Each sample is declared as an
# output, so a missing file fails the build.
rule megadepth:
    input:
        alignment=lambda wildcards: alignment(wildcards.sample),
        index=lambda wildcards: index(wildcards.sample)
    output:
        "data/{sample}.bw"
    threads: int(config["threads_per_job"])
    params:
        reference=reference_option("--fasta"),
        prefix="data/{sample}"
    shell:
        """
        megadepth {input.alignment:q} {params.reference} --bigwig --filter-out 1796 \
            --threads {threads} --prefix {params.prefix}
        mv {params.prefix}.all.bw {output:q}
        """


# The scans below run the rules bigwig_stat and bigwig_over10 of
# workflows/create/analysis.smk through the same timed_trials.py, with the same
# scripts, used in place. What
# differs: the data is in the working directory, threads are CORES, the cache is
# always cold, and trials come from the config.
#
# Everything a scan writes (per-sample tables, logs) goes in the working directory
# on local disk. Only the timing table goes to
# {results}/scans/bigwig.{analysis}.tsv, written by timed_trials.py after the last
# trial, so no Lustre write lands inside a timed scan. Scans run one at a time
# (resources: bench=1), because overlapping jobs evict and refill each other's
# page cache.

# One pyBigWig scan per sample, one command per line. The full scan reads whole
# chromosomes from each file's header, so it gets no BED. The commands go in a
# file because 1024 of them as arguments would pass the 128 KB limit.
rule bigwig_scan_commands:
    input:
        regions=lambda wc: [] if wc.analysis == "full_scan" else [REGIONS_BED]
    output:
        "scan/bigwig.{analysis}.commands"
    run:
        regions = " ".join(input.regions)
        with open(output[0], "w") as out:
            for sample in SAMPLES:
                out.write(f"{PYTHON} {BIGWIG_SCAN} data/{sample}.bw {regions}"
                          f" > scan/bigwig.{wildcards.analysis}/{sample}.tsv\n")


# xargs runs the command lines, CORES at a time, each with sh -c. The cache is
# controlled for the whole data/ folder, which holds only the .bw files.
rule bigwig_stat:
    input:
        bigwigs=BIGWIGS,
        commands="scan/bigwig.{analysis}.commands",
        scanner=BIGWIG_SCAN,
        timer=TIMED_TRIALS
    output:
        tsv=directory("scan/bigwig.{analysis}"),
        bench=RESULTS + "/scans/bigwig.{analysis}.tsv"
    threads: CORES
    resources:
        bench=1
    params:
        python=PYTHON,
        trials=TRIALS
    log:
        "logs/bigwig.{analysis}.log"
    shell:
        r"""
        mkdir -p {output.tsv:q}
        {params.python:q} {input.timer:q} --data data --cache cold \
            --trials {params.trials} --output {output.bench:q} -- \
            xargs -a {input.commands:q} -d '\n' -n 1 -P {threads} sh -c 2> {log:q}
        """


# Per-base fraction of samples with depth >= 10x (gnomAD's over_10), written as a
# bigWig by a plain script that uses a multiprocessing pool. Each worker reads one
# 1 Mb window from every sample's bigWig in turn, so its memory does not grow
# with N. The sample paths go on the command line (about 20 KB at N = 1024).
# --before removes the previous trial's output outside the timed window.
rule bigwig_over10:
    input:
        bigwigs=BIGWIGS,
        script=OVER10_BIGWIG,
        timer=TIMED_TRIALS
    output:
        out="scan/bigwig.over10.bw",
        bench=RESULTS + "/scans/bigwig.over10.tsv"
    threads: CORES
    resources:
        bench=1
    params:
        python=PYTHON,
        trials=TRIALS
    log:
        "logs/bigwig.over10.log"
    shell:
        """
        {params.python:q} {input.timer:q} --data data --cache cold \
            --trials {params.trials} --output {output.bench:q} --before 'rm -rf {output.out}' -- \
            {params.python:q} {input.script:q} {output.out:q} {threads} {input.bigwigs:q} \
            2> {log:q}
        """
