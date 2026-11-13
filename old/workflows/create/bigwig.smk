include: "samples.smk"

SELECTED = {
    s["sample_id"]: s for s in selected_samples(config["n_samples"])
}


rule all:
    input:
        expand("data/{sample}.bw", sample=SELECTED)


rule megadepth:
    input:
        alignment=lambda wildcards: SELECTED[wildcards.sample]["alignment"],
        index=lambda wildcards: SELECTED[wildcards.sample]["index"],
        reference=REFERENCE,
        reference_index=REFERENCE_INDEX
    output:
        "data/{sample}.bw"
    threads: 4
    params:
        prefix="data/{sample}"
    benchmark:
        "benchmarks/megadepth/{sample}.tsv"
    shell:
        """
        megadepth {input.alignment:q} --fasta {input.reference:q} \
            --bigwig --filter-out 1796 --threads {threads} --prefix {params.prefix:q}
        mv {params.prefix:q}.all.bw {output:q}
        """
