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
        {params.pbz:q} import bam --file-list {input.sources:q} --reference {params.reference:q} \
            --threads {threads} --track depth -c sample --codecs @{input.codecs:q} \
            --no-progress -o {output:q}
        """
