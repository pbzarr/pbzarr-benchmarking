# Cohort workflow for BigWig: one Megadepth run per sample, then the scans.
#   snakemake -s workflows/cohort/bigwig.smk build   (timed by the driver, with timed_trials.py)
#   snakemake -s workflows/cohort/bigwig.smk scans   (each scan timed by timed_trials.py)
# build: --cores JOBS * THREADS. scans: --cores <cores> --resources bench=1.
# See common.smk.

include: "common.smk"

wildcard_constraints:
    sample="[^/]+"

BIGWIG_SCAN = str(CREATE / "bigwig_scan.py")


rule build:
    input:
        BIGWIGS


rule scans:
    input:
        expand(RESULTS + "/scans/bigwig.{analysis}.tsv", analysis=["full_scan", "region_scan"])


# Each job declares threads_per_job threads, so with --cores JOBS * threads_per_job
# Snakemake runs JOBS of them at once. Each sample is declared as an output, so a
# missing file fails the build.
rule megadepth:
    input:
        cram=lambda wildcards: cram(wildcards.sample),
        crai=lambda wildcards: cram(wildcards.sample) + ".crai"
    output:
        "data/{sample}.bw"
    threads: int(config["threads_per_job"])
    params:
        reference=REFERENCE,
        prefix="data/{sample}"
    shell:
        """
        megadepth {input.cram:q} --fasta {params.reference:q} --bigwig --filter-out 1796 \
            --threads {threads} --prefix {params.prefix}
        mv {params.prefix}.all.bw {output:q}
        """


# The scans below run the rule bigwig_stat of workflows/create/analysis.smk
# through the same timed_trials.py, with the same script, used in place. What
# differs: the data is in the working directory, threads are CORES, the cache is
# always cold, and trials come from the config.
#
# Everything a scan writes (per-sample tables, logs) goes in the working directory
# on local disk. Only the timing table goes to
# {results}/scans/bigwig.{analysis}.tsv, written by timed_trials.py after the last
# trial, so no Lustre write lands inside a timed scan. Scans run one at a time
# (resources: bench=1), because overlapping jobs evict and refill each other's
# page cache.

# One pyBigWig scan per sample, run in parallel with CORES processes in total. The
# cache is controlled for the whole data/ folder, which holds only the .bw files.
rule bigwig_stat:
    input:
        bigwigs=BIGWIGS,
        scanner=BIGWIG_SCAN,
        regions=lambda wc: [] if wc.analysis == "full_scan" else [REGIONS_BED],
        timer=TIMED_TRIALS
    output:
        tsv=directory("scan/bigwig.{analysis}"),
        bench=RESULTS + "/scans/bigwig.{analysis}.tsv"
    threads: CORES
    resources:
        bench=1
    params:
        python=PYTHON,
        trials=TRIALS,
        regions=lambda wc: "/dev/null" if wc.analysis == "full_scan" else REGIONS_BED
    log:
        "logs/bigwig.{analysis}.log"
    shell:
        r"""
        test {threads} -eq {CORES}
        mkdir -p {output.tsv:q}
        scan=$(cat <<'EOF'
printf '%s\0' {input.bigwigs:q} |
    xargs -0 -n 1 -P {threads} sh -c '
        python="$1"
        scanner="$2"
        analysis="$3"
        regions="$4"
        prefix="$5/$(basename "$6")"
        if [ "$analysis" = full_scan ]; then
            "$python" "$scanner" "$6"
        else
            "$python" "$scanner" "$6" "$regions"
        fi > "$prefix.tsv"
    ' sh {params.python:q} {input.scanner:q} {wildcards.analysis:q} {params.regions:q} \
        {output.tsv:q}
EOF
)
        {params.python:q} {input.timer:q} --data data --cache cold \
            --trials {params.trials} --output {output.bench:q} -- \
            sh -c "$scan" 2> {log:q}
        """
