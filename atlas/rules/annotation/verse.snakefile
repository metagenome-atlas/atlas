rule counts_per_region:
    input:
        gtf = "{sample}/%s/annotation/orfs/{sample}.gtf" % ASSEMBLER,
        bam = "{sample}/%s/annotation/{sample}.bam" % ASSEMBLER,
        bai = "{sample}/%s/annotation/{sample}.bam.bai" % ASSEMBLER
    output:
        summary = "{sample}/%s/annotation/orfs/{sample}.CDS.summary.txt" % ASSEMBLER,
        counts = "{sample}/%s/annotation/orfs/{sample}.CDS.txt" % ASSEMBLER
    params:
        min_read_overlap = config["annotation"].get("minimum_overlap", 20)
    log:
        "{sample}/%s/logs/counts_per_region.log" % ASSEMBLER
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} verse -T {threads} --minReadOverlap {params.min_read_overlap} \
               --singleEnd -t CDS -z 1 -a {input.gtf} \
               -o {wildcards.sample}/{ASSEMBLER}/annotation/orfs/{wildcards.sample} \
               {input.bam} > {log}"""
