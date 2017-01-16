rule spades:
    input:
        "{sample}/quality_control/{NORMALIZATION}/{sample}_pe.fastq.gz"
    output:
        temp("{sample}/{ASSEMBLER}/contigs.fasta")
    params:
        # memory = config["assembly"].get("memory", 0.90)
        k = config["assembly"].get("spades_k", "auto"),
        outdir = lambda wc: "{sample}/{assembler}".format(sample=wc.sample, assembler=ASSEMBLER)
    log:
        "{sample}/{ASSEMBLER}/spades.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} spades.py -t {threads} -o {params.outdir} --meta --12 {input}"""


rule rename_spades_output:
    input:
        "{sample}/{ASSEMBLER}/contigs.fasta"
    output:
        "{sample}/{ASSEMBLER}/{sample}_prefilter_contigs.fasta"
    shell:
        "{SHPFXS} cp {input} {output}"
