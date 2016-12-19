rule spades:
    input:
        "results/{eid}/{sample}/quality_control/%s/{sample}_pe.fastq.gz" % NORMALIZATION
    output:
        temp("results/{eid}/{sample}/%s/contigs.fasta" % ASSEMBLER)
    params:
        # memory = config["assembly"].get("memory", 0.90)
        k = config["assembly"].get("spades_k", "auto"),
        outdir = lambda wc: "results/%s/%s/%s" % (wc.eid, wc.sample, ASSEMBLER)
    log:
        "results/{eid}/{sample}/%s/spades.log" % ASSEMBLER
    shadow:
        "shallow"
    threads:
        config.get("threads", 1)
    # shadow makes threads a little bit different
    shell:
        """{SHPFXM} spades.py -t {config[threads]} -o {params.outdir} --meta --12 {input}"""


rule rename_spades_output:
    input:
        "results/{eid}/{sample}/%s/contigs.fasta" % ASSEMBLER
    output:
        "results/{eid}/{sample}/%s/{sample}_prefilter_contigs.fasta" % ASSEMBLER
    shell:
        "{SHPFXS} cp {input} {output}"
