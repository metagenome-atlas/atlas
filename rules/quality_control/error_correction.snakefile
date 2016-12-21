rule error_correction:
    input:
        "{sample}/quality_control/quality_filter/{sample}_pe.fastq.gz"
    output:
        "{sample}/quality_control/error_correction/{sample}_pe.fastq.gz"
    log:
        "{sample}/logs/{sample}_error_correction.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} tadpole.sh in={input} out={output} mode=correct threads={threads} \
               ecc=t ecco=t 2> {log}"""
