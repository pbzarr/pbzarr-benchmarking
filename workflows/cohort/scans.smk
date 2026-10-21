# The scans over a finished cohort build: full scan, region scan, and over10
# (PBZ and D4 only). Included by every cohort workflow; each workflow asks only
# for its own tool's tables through its `scans` target.
#
# The commands are the rules pbz_stat, d4_stat, bigwig_stat, pbz_over10 and
# d4_over10 of workflows/create/analysis.smk, run through the same
# timed_trials.py with the same scripts, used in place. What differs: the data
# is in the working directory, threads are CORES, the cache is always cold, and
# trials come from the config.
#
# Everything a scan writes (stat tables, over10 outputs, BigWig per-sample
# tables, logs) goes in the working directory on local disk. Only the timing
# table goes to {results}/scans/{tool}.{analysis}.tsv, written by timed_trials.py
# after the last trial, so no Lustre write lands inside a timed scan.
#
# Scans run one at a time (resources: bench=1), because overlapping jobs evict
# and refill each other's page cache.

CREATE = ROOT / "workflows/create"
TIMED_TRIALS = str(CREATE / "timed_trials.py")
BIGWIG_SCAN = str(CREATE / "bigwig_scan.py")
OVER10_PBZ = str(CREATE / "over10_pbz.py")
OVER10_D4 = str(CREATE / "over10_d4.py")
PROPOSED_API = CREATE / "proposed_api"
REGIONS_BED = CFG["regions_bed"]
PYTHON = sys.executable
TRIALS = int(CFG["scan_trials"])

wildcard_constraints:
    analysis="full_scan|region_scan"


rule full_scan_bed:
    input:
        CFG["chrom_sizes"]
    output:
        "full_scan.bed"
    shell:
        r"""awk 'BEGIN {{OFS="\t"}} $2 > 0 {{print $1, 0, $2}}' {input:q} > {output:q}"""


rule pbz_stat:
    input:
        store=PBZ_STORE,
        regions=lambda wc: "full_scan.bed" if wc.analysis == "full_scan" else REGIONS_BED,
        timer=TIMED_TRIALS
    output:
        tsv="scan/pbz.{analysis}.tsv",
        bench=RESULTS + "/scans/pbz.{analysis}.tsv"
    threads: CORES
    resources:
        bench=1
    params:
        pbz=CFG["pbz"],
        python=PYTHON,
        trials=TRIALS
    log:
        "logs/pbz.{analysis}.log"
    shell:
        """
        {params.python:q} {input.timer:q} --data {input.store:q} --cache cold \
            --trials {params.trials} --output {output.bench:q} -- \
            {params.pbz:q} stat {input.store:q} depth \
            -s mean -t {threads} -r {input.regions:q} -o {output.tsv:q} 2> {log:q}
        """


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
        r"""
        scan=$(cat <<'EOF'
d4tools stat {params.scan_flag} -s mean -H -t {threads} \
    -r {input.regions:q} {input.d4:q} > {output.tsv:q}
EOF
)
        {params.python:q} {input.timer:q} --data {input.d4:q} --cache cold \
            --trials {params.trials} --output {output.bench:q} -- \
            sh -c "$scan" 2> {log:q}
        """


# Per-base fraction of samples with depth >= 10x (gnomAD's over_10), written in
# each format's own file type by one plain script per format: a dask cluster
# for PBZ (pbzarr is dask-backed), a multiprocessing pool for D4.
# --before removes the previous trial's output outside the timed window.

rule pbz_over10:
    input:
        store=PBZ_STORE,
        script=OVER10_PBZ,
        shim=str(PROPOSED_API / "pbzarr/__init__.py"),
        cluster=str(PROPOSED_API / "pbzarr/cluster.py"),
        timer=TIMED_TRIALS
    output:
        out=directory("scan/pbz.over10.pbz"),
        bench=RESULTS + "/scans/pbz.over10.tsv"
    threads: CORES
    resources:
        bench=1
    params:
        python=PYTHON,
        trials=TRIALS,
        proposed_api=str(PROPOSED_API)
    log:
        "logs/pbz.over10.log"
    shell:
        """
        # over10_pbz.py imports the proposed-API shim as pbzarr
        PYTHONPATH={params.proposed_api:q} \
        {params.python:q} {input.timer:q} --data {input.store:q} --cache cold \
            --trials {params.trials} --output {output.bench:q} --before 'rm -rf {output.out}' -- \
            {params.python:q} {input.script:q} {input.store:q} {output.out:q} {threads} \
            2> {log:q}
        """


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
        {params.python:q} {input.timer:q} --data {input.d4:q} --cache cold \
            --trials {params.trials} --output {output.bench:q} --before 'rm -rf {output.out}' -- \
            {params.python:q} {input.script:q} {input.d4:q} {output.out:q} {threads} \
            2> {log:q}
        """
