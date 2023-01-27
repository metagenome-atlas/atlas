bin_quality_input_folder = "{sample}/binning/{binner}/bins"


### Checkm2 ###


rule checkm2_download_db:
    output:
        directory(f"{DBDIR}/CheckM2"),
    conda:
        "../envs/checkm2.yaml"
    threads: 1
    log:
        "logs/download/checkm2.log",
    resources:
        time=config["runtime"]["long"],
    shell:
        " checkm2 database --download --path {output} "
        " &>> {log}"


rule run_checkm2:
    input:
        fasta_dir=bin_quality_input_folder,
        db=rules.checkm2_download_db.output,
    output:
        directory("{sample}/binning/{binner}/bin_quality/checkm2"),
    conda:
        "../envs/checkm2.yaml"
    threads: config["threads"]
    log:
        "{sample}/logs/binning/{binner}/checkm2.log",
    benchmark:
        "logs/benchmarks/checkm2/{sample}_{binner}.tsv"
    resources:
        time=int(config["runtime"]["default"]),
        mem_mb=config["mem"],
    shell:
        " checkm2 predict "
        " --threads {threads} "
        " --force "
        " -x .fasta "
        " --input {input.fasta_dir} "
        " --output-directory {output} "
        " &> {log} "

localrules:
    build_bin_report,
    combine_bin_stats,


rule combine_bin_stats:
    input:
        completeness_files=expand(
            "{sample}/binning/{{binner}}/bin_quality/checkm2", sample=SAMPLES
        ),
    output:
        bin_table="reports/genomic_bins_{binner}.tsv",
    params:
        samples=SAMPLES,
    log:
        "logs/binning/combine_stats_{binner}.log",
    script:
        "../scripts/combine_checkm2.py"


rule build_bin_report:
    input:
        bin_table="reports/genomic_bins_{binner}.tsv",
    output:
        report="reports/bin_report_{binner}.html",
    conda:
        "../envs/report.yaml"
    log:
        "logs/binning/report_{binner}.log",
    script:
        "../report/bin_report.py"
