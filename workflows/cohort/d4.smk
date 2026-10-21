# Cohort workflow for D4: one mosdepth run per sample, then one `d4tools merge`,
# then the scans.
#   snakemake -s workflows/cohort/d4.smk build   (timed by the driver's Slurm step)
#   snakemake -s workflows/cohort/d4.smk scans   (each scan timed by timed_trials.py)

include: "common.smk"
include: "scans.smk"

wildcard_constraints:
    sample="[^/]+"

EVICT = str(ROOT / "workflows/lib/evict.py")
PER_SAMPLE = [f"tmp/{sample}.per-base.d4" for sample in SAMPLES]


rule build:
    input:
        D4_FILE


rule scans:
    input:
        expand(RESULTS + "/scans/d4.{analysis}.tsv", analysis=["full_scan", "region_scan", "over10"])


# Up to `pack` jobs run at once, each with threads_per_job threads. Each sample
# is declared as an output, so a missing file fails the build.
rule mosdepth:
    input:
        cram=lambda wildcards: cram(wildcards.sample),
        crai=lambda wildcards: cram(wildcards.sample) + ".crai"
    output:
        temp("tmp/{sample}.per-base.d4")
    threads: int(config["threads_per_job"])
    resources:
        pack=1
    params:
        reference=REFERENCE,
        prefix="tmp/{sample}"
    shell:
        "mosdepth -t {threads} --fasta {params.reference:q} --d4 {params.prefix} {input.cram:q}"


# Runs between the last mosdepth job and the merge, so the merge reads its inputs
# cold. The driver reads d4_per_sample_bytes.txt for the peak disk use.
rule evict_per_sample:
    input:
        PER_SAMPLE
    output:
        "tmp/.evicted"
    params:
        evict=EVICT,
        results=RESULTS
    shell:
        """
        du -cb {input:q} | tail -1 | cut -f1 > {params.results:q}/d4_per_sample_bytes.txt
        sync
        python {params.evict:q} {input:q}
        touch {output:q}
        """


# The only step these workflows time themselves. The row goes to merge.tsv in
# the results folder after the timed command has finished.
rule merge:
    input:
        per_sample=PER_SAMPLE,
        evicted="tmp/.evicted"
    output:
        D4_FILE
    params:
        sources=[f"tmp/{sample}.per-base.d4:{sample}" for sample in SAMPLES],
        results=RESULTS,
        n=config["n_samples"],
        build=config["build"]
    shell:
        TIMED("d4tools merge {params.sources} {output}", "merge.time") + r"""
        read wall user sys rss < merge.time
        bytes=$(du -sb {output:q} | cut -f1)
        table={params.results:q}/merge.tsv
        [ -e "$table" ] || printf 'n\tbuild\twall_s\tuser_s\tsys_s\tmax_rss_kb\tout_bytes\n' > "$table"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' {params.n} {params.build} "$wall" "$user" "$sys" "$rss" "$bytes" >> "$table"
        """
