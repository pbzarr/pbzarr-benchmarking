# Cohort workflow for PBZ: one `pbz import` over all samples, then the scans.
#   snakemake -s workflows/cohort/pbz.smk build   (timed by the driver's Slurm step)
#   snakemake -s workflows/cohort/pbz.smk scans   (each scan timed by timed_trials.py)

include: "common.smk"
include: "scans.smk"


rule build:
    input:
        PBZ_STORE


rule scans:
    input:
        expand(RESULTS + "/scans/pbz.{analysis}.tsv", analysis=["full_scan", "region_scan", "over10"])


# One process imports every sample, in manifest order, at CORES threads.
# -c sample is needed even for a single sample.
rule import:
    input:
        crams=[cram(sample) for sample in SAMPLES],
        crais=[cram(sample) + ".crai" for sample in SAMPLES],
        codecs=CFG["codecs"]
    output:
        directory(PBZ_STORE)
    threads: CORES
    params:
        pbz=CFG["pbz"],
        reference=REFERENCE,
        sources=[f"{cram(sample)}:{sample}" for sample in SAMPLES]
    shell:
        """
        {params.pbz:q} import bam {params.sources:q} --reference {params.reference:q} \
            --threads {threads} --track depth -c sample --codecs @{input.codecs:q} \
            --no-progress -o {output:q}
        """
