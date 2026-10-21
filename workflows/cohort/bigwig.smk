# Cohort workflow for BigWig: one Megadepth run per sample, then the scans.
#   snakemake -s workflows/cohort/bigwig.smk build   (timed by the driver's Slurm step)
#   snakemake -s workflows/cohort/bigwig.smk scans   (each scan timed by timed_trials.py)

include: "common.smk"
include: "scans.smk"

wildcard_constraints:
    sample="[^/]+"


rule build:
    input:
        BIGWIGS


rule scans:
    input:
        expand(RESULTS + "/scans/bigwig.{analysis}.tsv", analysis=["full_scan", "region_scan"])


# Up to `pack` jobs run at once, each with threads_per_job threads. Each sample
# is declared as an output, so a missing file fails the build.
rule megadepth:
    input:
        cram=lambda wildcards: cram(wildcards.sample),
        crai=lambda wildcards: cram(wildcards.sample) + ".crai"
    output:
        "data/{sample}.bw"
    threads: int(config["threads_per_job"])
    resources:
        pack=1
    params:
        reference=REFERENCE,
        prefix="data/{sample}"
    shell:
        """
        megadepth {input.cram:q} --fasta {params.reference:q} --bigwig --filter-out 1796 \
            --threads {threads} --prefix {params.prefix}
        mv {params.prefix}.all.bw {output:q}
        """
