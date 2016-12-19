rule fastqc:
    input:
        "results/{eid}/{sample}/quality_control/error_correction/{sample}_pe.fastq.gz"
    output:
        "results/{eid}/{sample}/quality_control/fastqc/{sample}_pe_fastqc.zip",
        "results/{eid}/{sample}/quality_control/fastqc/{sample}_pe_fastqc.html"
    params:
        output_dir = lambda wc: "results/%s/%s/quality_control/fastqc/" % (wc.eid, wc.sample)
    threads:
        config.get("threads", 1)
    shell:
        "{SHPFXM} fastqc -t {threads} -f fastq -o {params.output_dir} {input}"
