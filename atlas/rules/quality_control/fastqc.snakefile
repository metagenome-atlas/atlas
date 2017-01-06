rule fastqc:
    input:
        "{sample}/quality_control/error_correction/{sample}_pe.fastq.gz"
    output:
        "{sample}/quality_control/fastqc/{sample}_pe_fastqc.zip",
        "{sample}/quality_control/fastqc/{sample}_pe_fastqc.html"
    params:
        output_dir = lambda wc: "{sample}/quality_control/fastqc/".format(sample=wc.sample)
    threads:
        config.get("threads", 1)
    shell:
        "{SHPFXM} fastqc -t {threads} -f fastq -o {params.output_dir} {input}"
