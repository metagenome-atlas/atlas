rule megahit:
    input:
        "{sample}/quality_control/%s/{sample}_pe.fastq.gz" % NORMALIZATION
    output:
        temp("{sample}/%s/{sample}_prefilter.contigs.fa" % ASSEMBLER)
    params:
        memory = config["assembly"].get("memory", 0.90),
        min_count = config["assembly"].get("minimum_count", 2),
        k_min = config["assembly"].get("kmer_min", 21),
        k_max = config["assembly"].get("kmer_max", 121),
        k_step = config["assembly"].get("kmer_step", 20),
        merge_level = config["assembly"].get("merge_level", "20,0.98"),
        prune_level = config["assembly"].get("prune_level", 2),
        low_local_ratio = config["assembly"].get("low_local_ratio", 0.2),
        min_contig_len = config["assembly"].get("minimum_contig_length", 200),
        outdir = lambda wc: "{sample}/{assembler}".format(sample=wc.sample, assembler=ASSEMBLER)
    log:
        "{sample}/%s/{sample}.log" % ASSEMBLER
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} megahit --num-cpu-threads {threads} --12 {input} --continue \
               --k-min {params.k_min} --k-max {params.k_max} --k-step {params.k_step} \
               --out-dir {params.outdir} --out-prefix {wildcards.sample}_prefilter \
               --min-contig-len {params.min_contig_len} --min-count {params.min_count} \
               --merge-level {params.merge_level} --prune-level {params.prune_level} \
               --low-local-ratio {params.low_local_ratio}"""


rule rename_megahit_output:
    input:
        "{sample}/%s/{sample}_prefilter.contigs.fa" % ASSEMBLER
    output:
        "{sample}/%s/{sample}_prefilter_contigs.fasta" % ASSEMBLER
    shell:
        "{SHPFXS} cp {input} {output}"
