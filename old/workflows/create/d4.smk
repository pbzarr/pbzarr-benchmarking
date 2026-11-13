include: "samples.smk"

SELECTED = {
    s["sample_id"]: s for s in selected_samples(config["n_samples"])
}
D4_SPARSE = config.get("d4_sparse", "rna" in config["dataset"].lower())
D4_DEFLATE = config.get("d4_deflate", True)
if type(D4_SPARSE) is not bool or type(D4_DEFLATE) is not bool:
    raise ValueError("d4_sparse and d4_deflate must be YAML booleans (true or false)")


rule all:
    input:
        "data/cohort.d4"


if D4_SPARSE:
    rule encode_sparse_d4:
        input:
            alignment=lambda wildcards: SELECTED[wildcards.sample]["alignment"],
            index=lambda wildcards: SELECTED[wildcards.sample]["index"],
            reference=REFERENCE,
            reference_index=REFERENCE_INDEX,
            dictionary=str(WORKFLOW_DIR / "d4-zero.dict")
        output:
            temp("tmp/d4/{sample}.d4")
        threads: 4
        params:
            compression="-z" if D4_DEFLATE else ""
        benchmark:
            "benchmarks/encode_sparse_d4/{sample}.tsv"
        shell:
            "d4tools create -t {threads} {params.compression} "
            "--dict-file {input.dictionary:q} -r {input.reference:q} "
            "{input.alignment:q} {output:q}"
else:
    rule mosdepth:
        input:
            alignment=lambda wildcards: SELECTED[wildcards.sample]["alignment"],
            index=lambda wildcards: SELECTED[wildcards.sample]["index"],
            reference=REFERENCE,
            reference_index=REFERENCE_INDEX
        output:
            temp(directory("tmp/{sample}"))
        threads: 4
        benchmark:
            "benchmarks/mosdepth/{sample}.tsv"
        shell:
            "mkdir -p {output:q} && "
            "mosdepth -t 4 --fasta {input.reference:q} --d4 "
            "{output:q}/coverage {input.alignment:q}"


rule merge:
    input:
        expand("tmp/d4/{sample}.d4", sample=SELECTED) if D4_SPARSE else
        expand("tmp/{sample}", sample=SELECTED)
    output:
        "data/cohort.d4"
    params:
        sources=[f"tmp/d4/{sample}.d4:{sample}" for sample in SELECTED] if D4_SPARSE else [
            f"tmp/{sample}/coverage.per-base.d4:{sample}"
            for sample in SELECTED
        ]
    benchmark:
        "benchmarks/merge.tsv"
    shell:
        "d4tools merge {params.sources:q} {output:q}"
