import sys


PBZ_BINARY = "/uufs/chpc.utah.edu/common/HIPAA/u1364833/pbzarr/target/release/pbz"
REGIONS_BED = str(ROOT / config["regions_bed"])
FULL_SCAN_BED = RESULTS + "/full_scan.bed"
BIGWIG_SCAN = str(WORKFLOW_DIR / "bigwig_scan.py")
TIMED_TRIALS = str(WORKFLOW_DIR / "timed_trials.py")
# Back-to-back timed trials per scan job. cold evicts the scanned data from
# the page cache before every trial; warm runs one untimed pass first.
ANALYSIS_TRIALS = int(config.get("analysis_trials", 5))
# Every timed job takes the one scan_slot, so run with --resources scan_slot=1:
# overlapping jobs evict and refill each other's page cache.

wildcard_constraints:
    analysis="full_scan|region_scan",
    cache="cold|warm"


rule full_scan_regions:
    input:
        CHROM_SIZES
    output:
        FULL_SCAN_BED
    shell:
        r"""awk 'BEGIN {{OFS="\t"}} $2 > 0 {{print $1, 0, $2}}' {input:q} > {output:q}"""


rule pbz_stat:
    input:
        data=RESULTS + "/pbz/t{threads}/n{n_samples}/r1/data",
        regions=lambda wc: FULL_SCAN_BED if wc.analysis == "full_scan" else REGIONS_BED,
        timer=TIMED_TRIALS
    output:
        tsv=RESULTS + "/pbz/t{threads}/n{n_samples}/r{replicate}/{analysis}.{cache}.tsv",
        bench=RESULTS + "/pbz/t{threads}/n{n_samples}/r{replicate}/benchmarks/{analysis}.{cache}.tsv"
    threads: lambda wc: int(wc.threads)
    resources:
        scan_slot=1
    params:
        pbz_binary=PBZ_BINARY,
        python=sys.executable,
        trials=ANALYSIS_TRIALS
    log:
        RESULTS + "/pbz/t{threads}/n{n_samples}/r{replicate}/{analysis}.{cache}.log"
    shell:
        """
        test {threads} -eq {wildcards.threads}
        {params.python:q} {input.timer:q} --data {input.data:q} --cache {wildcards.cache} \
            --trials {params.trials} --output {output.bench:q} -- \
            {params.pbz_binary:q} stat {input.data:q}/cohort.pbz depth \
            -s mean -t {threads} -r {input.regions:q} -o {output.tsv:q} 2> {log:q}
        """


rule bigwig_stat:
    input:
        data=RESULTS + "/bigwig/t{threads}/n{n_samples}/r1/data",
        scanner=BIGWIG_SCAN,
        regions=lambda wc: [] if wc.analysis == "full_scan" else [REGIONS_BED],
        timer=TIMED_TRIALS
    output:
        tsv=directory(RESULTS + "/bigwig/t{threads}/n{n_samples}/r{replicate}/{analysis}.{cache}"),
        bench=RESULTS + "/bigwig/t{threads}/n{n_samples}/r{replicate}/benchmarks/{analysis}.{cache}.tsv"
    threads: lambda wc: int(wc.threads)
    resources:
        scan_slot=1
    params:
        python=sys.executable,
        trials=ANALYSIS_TRIALS,
        regions=lambda wc: "/dev/null" if wc.analysis == "full_scan" else REGIONS_BED,
        bigwigs=lambda wc, input: [
            f"{input.data}/{s['sample_id']}.bw" for s in selected_samples(wc.n_samples)
        ]
    log:
        RESULTS + "/bigwig/t{threads}/n{n_samples}/r{replicate}/{analysis}.{cache}.log"
    shell:
        r"""
        test {threads} -eq {wildcards.threads}
        mkdir -p {output.tsv:q}
        scan=$(cat <<'EOF'
printf '%s\0' {params.bigwigs:q} |
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
        {params.python:q} {input.timer:q} --data {input.data:q} --cache {wildcards.cache} \
            --trials {params.trials} --output {output.bench:q} -- \
            sh -c "$scan" 2> {log:q}
        """


rule d4_stat:
    input:
        data=RESULTS + "/d4/t{threads}/n{n_samples}/r1/data",
        regions=lambda wc: FULL_SCAN_BED if wc.analysis == "full_scan" else REGIONS_BED,
        timer=TIMED_TRIALS
    output:
        tsv=RESULTS + "/d4/t{threads}/n{n_samples}/r{replicate}/{analysis}.{cache}.tsv",
        bench=RESULTS + "/d4/t{threads}/n{n_samples}/r{replicate}/benchmarks/{analysis}.{cache}.tsv"
    threads: lambda wc: int(wc.threads)
    resources:
        scan_slot=1
    params:
        python=sys.executable,
        trials=ANALYSIS_TRIALS,
        scan_flag=lambda wc: "--no-index" if wc.analysis == "full_scan" else ""
    log:
        RESULTS + "/d4/t{threads}/n{n_samples}/r{replicate}/{analysis}.{cache}.log"
    shell:
        r"""
        test {threads} -eq {wildcards.threads}
        scan=$(cat <<'EOF'
d4tools stat {params.scan_flag} -s mean -H -t {threads} \
    -r {input.regions:q} {input.data:q}/cohort.d4 > {output.tsv:q}
EOF
)
        {params.python:q} {input.timer:q} --data {input.data:q} --cache {wildcards.cache} \
            --trials {params.trials} --output {output.bench:q} -- \
            sh -c "$scan" 2> {log:q}
        """


# Per-base fraction of samples with depth >= 10x (gnomAD's over_10), written in
# each format's own file type by one plain script per format: a dask cluster
# for PBZ (pbzarr is dask-backed), a multiprocessing pool for BigWig and D4.
# --before removes the previous trial's output outside the timed window.
OVER10 = {tool: str(WORKFLOW_DIR / f"over10_{tool}.py") for tool in ("pbz", "d4", "bigwig")}


rule pbz_over10:
    input:
        data=RESULTS + "/pbz/t{threads}/n{n_samples}/r1/data",
        script=OVER10["pbz"],
        shim=str(WORKFLOW_DIR / "proposed_api/pbzarr/__init__.py"),
        cluster=str(WORKFLOW_DIR / "proposed_api/pbzarr/cluster.py"),
        timer=TIMED_TRIALS
    output:
        out=directory(RESULTS + "/pbz/t{threads}/n{n_samples}/r{replicate}/over10.{cache}.pbz"),
        bench=RESULTS + "/pbz/t{threads}/n{n_samples}/r{replicate}/benchmarks/over10.{cache}.tsv"
    threads: lambda wc: int(wc.threads)
    resources:
        scan_slot=1
    params:
        python=sys.executable,
        trials=ANALYSIS_TRIALS,
        proposed_api=str(WORKFLOW_DIR / "proposed_api")
    log:
        RESULTS + "/pbz/t{threads}/n{n_samples}/r{replicate}/over10.{cache}.log"
    shell:
        """
        test {threads} -eq {wildcards.threads}
        # over10_pbz.py imports the proposed-API shim as pbzarr. The malloc tunables let
        # glibc reuse freed chunk buffers instead of mmapping fresh pages for each one.
        PYTHONPATH={params.proposed_api:q} \
        GLIBC_TUNABLES=glibc.malloc.mmap_threshold=33554432:glibc.malloc.trim_threshold=17179869184 \
        {params.python:q} {input.timer:q} --data {input.data:q} --cache {wildcards.cache} \
            --trials {params.trials} --output {output.bench:q} --before 'rm -rf {output.out}' -- \
            {params.python:q} {input.script:q} {input.data:q}/cohort.pbz {output.out:q} {threads} \
            2> {log:q}
        """


rule d4_over10:
    input:
        data=RESULTS + "/d4/t{threads}/n{n_samples}/r1/data",
        script=OVER10["d4"],
        timer=TIMED_TRIALS
    output:
        out=RESULTS + "/d4/t{threads}/n{n_samples}/r{replicate}/over10.{cache}.d4",
        bench=RESULTS + "/d4/t{threads}/n{n_samples}/r{replicate}/benchmarks/over10.{cache}.tsv"
    threads: lambda wc: int(wc.threads)
    resources:
        scan_slot=1
    params:
        python=sys.executable,
        trials=ANALYSIS_TRIALS
    log:
        RESULTS + "/d4/t{threads}/n{n_samples}/r{replicate}/over10.{cache}.log"
    shell:
        """
        test {threads} -eq {wildcards.threads}
        {params.python:q} {input.timer:q} --data {input.data:q} --cache {wildcards.cache} \
            --trials {params.trials} --output {output.bench:q} --before 'rm -rf {output.out}' -- \
            {params.python:q} {input.script:q} {input.data:q}/cohort.d4 {output.out:q} {threads} \
            2> {log:q}
        """


rule bigwig_over10:
    input:
        data=RESULTS + "/bigwig/t{threads}/n{n_samples}/r1/data",
        script=OVER10["bigwig"],
        timer=TIMED_TRIALS
    output:
        out=RESULTS + "/bigwig/t{threads}/n{n_samples}/r{replicate}/over10.{cache}.bw",
        bench=RESULTS + "/bigwig/t{threads}/n{n_samples}/r{replicate}/benchmarks/over10.{cache}.tsv"
    threads: lambda wc: int(wc.threads)
    resources:
        scan_slot=1
    params:
        python=sys.executable,
        trials=ANALYSIS_TRIALS,
        bigwigs=lambda wc, input: [
            f"{input.data}/{s['sample_id']}.bw" for s in selected_samples(wc.n_samples)
        ]
    log:
        RESULTS + "/bigwig/t{threads}/n{n_samples}/r{replicate}/over10.{cache}.log"
    shell:
        """
        test {threads} -eq {wildcards.threads}
        {params.python:q} {input.timer:q} --data {input.data:q} --cache {wildcards.cache} \
            --trials {params.trials} --output {output.bench:q} --before 'rm -rf {output.out}' -- \
            {params.python:q} {input.script:q} {output.out:q} {threads} {params.bigwigs:q} \
            2> {log:q}
        """
