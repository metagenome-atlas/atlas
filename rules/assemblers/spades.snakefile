rule spades:
    input:
        "{sample}/quality_control/%s/{sample}_pe.fastq.gz" % NORMALIZATION
    output:
        temp("{sample}/%s/contigs.fasta" % ASSEMBLER)
    params:
        # memory = config["assembly"].get("memory", 0.90)
        k = config["assembly"].get("spades_k", "auto"),
        outdir = lambda wc: "{sample}/{assembler}".format(sample=wc.sample, assembler=ASSEMBLER)
    log:
        "{sample}/%s/spades.log" % ASSEMBLER
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} spades.py -t {threads} -o {params.outdir} --meta --12 {input}"""


rule rename_spades_output:
    input:
        "{sample}/%s/contigs.fasta" % ASSEMBLER
    output:
        "{sample}/%s/{sample}_prefilter_contigs.fasta" % ASSEMBLER
    shell:
        "{SHPFXS} cp {input} {output}"
