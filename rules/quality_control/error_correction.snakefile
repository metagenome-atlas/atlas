rule error_correction:
    input:
        "results/{eid}/{sample}/quality_control/quality_filter/{sample}_pe.fastq.gz"
    output:
        "results/{eid}/{sample}/quality_control/error_correction/{sample}_pe.fastq.gz"
    log:
        "results/{eid}/{sample}/logs/{sample}_error_correction.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} tadpole.sh in={input} out={output} mode=correct threads={threads} \
               ecc=t ecco=t 2> {log}"""
