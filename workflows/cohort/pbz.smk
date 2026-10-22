# Cohort workflow for PBZ: one `pbz import` over all samples, then the scans.
#   snakemake -s workflows/cohort/pbz.smk build   (timed by the driver, with timed_trials.py)
#   snakemake -s workflows/cohort/pbz.smk scans   (each scan timed by timed_trials.py)
# --cores is the config's `cores` for both; the scans also need --resources bench=1.
# See common.smk.

include: "common.smk"

PROPOSED_API = CREATE / "proposed_api"
OVER10_PBZ = str(CREATE / "over10_pbz.py")


rule build:
    input:
        PBZ_STORE


rule scans:
    input:
        expand(RESULTS + "/scans/pbz.{analysis}.tsv", analysis=["full_scan", "region_scan", "over10"])


# The sources go in a file because 1024 CRAM paths as arguments would nearly fill
# the 128 KB limit on one command-line argument. One line per sample, in
# manifest order: the CRAM path, a tab, the sample name. The file is outside
# data/ so it does not count towards the output size.
rule sources:
    output:
        "sources.tsv"
    run:
        with open(output[0], "w") as out:
            for sample in SAMPLES:
                out.write(f"{cram(sample)}\t{sample}\n")


# One process imports every sample, in file order, at CORES threads.
# -c sample is needed even for a single sample.
rule import:
    input:
        sources="sources.tsv",
        crams=[cram(sample) for sample in SAMPLES],
        crais=[cram(sample) + ".crai" for sample in SAMPLES],
        codecs=CFG["codecs"]
    output:
        directory(PBZ_STORE)
    threads: CORES
    params:
        pbz=CFG["pbz"],
        reference=REFERENCE
    shell:
        """
        test {threads} -eq {CORES}
        {params.pbz:q} import bam --file-list {input.sources:q} --reference {params.reference:q} \
            --threads {threads} --track depth -c sample --codecs @{input.codecs:q} \
            --no-progress -o {output:q}
        """


# The scans below run the rules pbz_stat and pbz_over10 of
# workflows/create/analysis.smk through the same timed_trials.py, with the same
# scripts, used in place. What differs: the data is in the working directory,
# threads are CORES, the cache is always cold, and trials come from the config.
#
# Everything a scan writes (stat tables, over10 outputs, logs) goes in the working
# directory on local disk. Only the timing table goes to
# {results}/scans/pbz.{analysis}.tsv, written by timed_trials.py after the last
# trial, so no Lustre write lands inside a timed scan. Scans run one at a time
# (resources: bench=1), because overlapping jobs evict and refill each other's
# page cache.

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
        test {threads} -eq {CORES}
        {params.python:q} {input.timer:q} --data {input.store:q} --cache cold \
            --trials {params.trials} --output {output.bench:q} -- \
            {params.pbz:q} stat {input.store:q} depth \
            -s mean -t {threads} -r {input.regions:q} -o {output.tsv:q} 2> {log:q}
        """


# Per-base fraction of samples with depth >= 10x (gnomAD's over_10), written as a
# PBZ store by a plain script that uses a dask cluster (pbzarr is dask-backed).
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
        test {threads} -eq {CORES}
        # over10_pbz.py imports the proposed-API shim as pbzarr
        PYTHONPATH={params.proposed_api:q} \
        {params.python:q} {input.timer:q} --data {input.store:q} --cache cold \
            --trials {params.trials} --output {output.bench:q} --before 'rm -rf {output.out}' -- \
            {params.python:q} {input.script:q} {input.store:q} {output.out:q} {threads} \
            2> {log:q}
        """
