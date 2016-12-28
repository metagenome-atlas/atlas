rule counts_per_region:
    input:
        gtf = "{sample}/annotation/orfs/{sample}.gtf",
        bam = "{sample}/annotation/{sample}.bam",
        bai = "{sample}/annotation/{sample}.bam.bai"
    output:
        summary = "{sample}/annotation/orfs/{sample}.CDS.summary.txt",
        counts = "{sample}/annotation/orfs/{sample}.CDS.txt"
    params:
        min_read_overlap = config["annotation"].get("minimum_overlap", 20)
    log:
        "{sample}/logs/counts_per_region.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} verse -T {threads} --minReadOverlap {params.min_read_overlap} \
               --singleEnd -t CDS -z 1 -a {input.gtf} \
               -o {wildcards.sample}/annotation/orfs/{wildcards.sample} \
               {input.bam} > {log}"""
