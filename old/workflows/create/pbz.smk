include: "samples.smk"
PBZ_BINARY = "/uufs/chpc.utah.edu/common/HIPAA/u1364833/pbzarr/target/release/pbz"
SELECTED = selected_samples(config["n_samples"])
CODECS = str(ROOT / config.get("codecs", "codecs.json"))  # 2D shards: 1,048,576 positions x 8 samples


rule all:
    input:
        "data/cohort.pbz"


rule import_pbz:
    input:
        alignments=[s["alignment"] for s in SELECTED],
        indexes=[s["index"] for s in SELECTED],
        reference=REFERENCE,
        reference_index=REFERENCE_INDEX,
        codecs=CODECS
    output:
        directory("data/cohort.pbz")
    threads: workflow.cores
    params:
        sources=[s["alignment"] + ":" + s["sample_id"] for s in SELECTED],
        pbz_binary = PBZ_BINARY
    shell:
        """
        {params.pbz_binary} import bam {params.sources:q} --reference {input.reference:q} \
            --in-flight 1000 --threads {threads} --codecs @{input.codecs:q} --track depth --no-progress -o {output:q}
        """
