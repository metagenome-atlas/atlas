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
        "{sample}/binning/{binner}/bin_quality/checkm2/quality_report.tsv",
        #faa=temp(directory("{sample}/binning/{binner}/bin_quality/checkm2/protein_files")),
    params:
        outdir= lambda wc, output: Path(output[0]).parent,
        lowmem = " --lowmem " if config["mem"] < 10 else ""
    conda:
        "../envs/checkm2.yaml"
    shadow:
        "minimal"
    threads: config["threads"]
    log:
        "{sample}/logs/binning/{binner}/checkm2.log",
        "{sample}/binning/{binner}/bin_quality/checkm2/checkm2.log"
    benchmark:
        "logs/benchmarks/checkm2/{sample}_{binner}.tsv"
    resources:
        time=int(config["runtime"]["default"]),
        mem_mb=config["mem"],
    shell:
        " checkm2 predict "
        " --threads {threads} "
        " {params.lowmem} "
        " --force "
        " --allmodels "
        " -x .fasta "
        " --tmpdir {resources.tmpdir} "
        " --input {input.fasta_dir} "
        " --output-directory {params.outdir} "
        " &> {log[0]} "


## GUNC ###


rule run_gunc:
    input:
        db=rules.download_gunc.output[0].format(**config),
        fasta_dir=bin_quality_input_folder,
    output:
        "{sample}/binning/{binner}/bin_quality/gunc.tsv",
    params:
        tmpdir=lambda wc: f"{config['tmpdir']}/gunc/{wc.sample}",
    conda:
        "../envs/gunc.yaml"
    threads: config["threads"]
    log:
        "{sample}/logs/binning/{binner}/gunc.log",
    benchmark:
        "logs/benchmarks/gunc/{sample}_{binner}.tsv"
    resources:
        time=int(config["runtime"]["default"]),
        mem_mb=config["mem"],
    shell:
        " mkdir -p {params.tmpdir}/ 2> {log} "
        " ; "
        " gunc run "
        " --threads {threads} "
        " --db_file {input.db} "
        " --input_dir {input.fasta_dir} "
        " --file_suffix .fasta "
        " --out_dir {params.tmpdir} &>> {log} "
        " ; "
        " mv {params.tmpdir}/*.tsv {output} 2>> {log}"


##### BUSCO  #########


rule run_busco:
    input:
        fasta_dir=bin_quality_input_folder,
        db=BUSCODIR,
    output:
        "{sample}/binning/{binner}/bin_quality/busco.tsv",
    params:
        tmpdir=lambda wc: f"{config['tmpdir']}/busco/{wc.sample}_{wc.binner}",
    conda:
        "../envs/busco.yaml"
    threads: config["threads"]
    log:
        "{sample}/logs/binning/{binner}/busco.log",
    benchmark:
        "logs/benchmarks/busco/{sample}_{binner}.tsv"
    resources:
        time=int(config["runtime"]["default"]),
        mem_mb=config["mem"],
    shell:
        " busco -i {input.fasta_dir} "
        " --auto-lineage-prok "
        " -m genome "
        " --out_path {params.tmpdir} "
        " -o output "
        " --download_path {input.db} "
        " -c {threads} "
        " --offline &> {log} "
        " ; "
        " mv {params.tmpdir}/output/batch_summary.txt {output} 2>> {log}"


# fetch also output/logs/busco.log


##### checkM  #########
rule run_checkm:
    input:
        touched_output="logs/checkm_init.txt",
        bins=bin_quality_input_folder,  # actualy path to fastas
    output:
        completeness="{sample}/binning/{binner}/bin_quality/checkm.tsv",
        taxonomy="{sample}/binning/{binner}/bin_quality/checkm_taxonomy.tsv",
    params:
        tmp_output_dir=lambda wc: f"{config['tmpdir']}/checkm/{wc.sample}_{wc.binner}",
        extension="fasta",
    conda:
        "../envs/checkm.yaml"
    threads: config["threads"]
    log:
        "{sample}/logs/binning/{binner}/checkm.log",
    benchmark:
        "logs/benchmarks/checkm_lineage_wf/{sample}_{binner}.tsv"
    resources:
        time=config["runtime"]["long"],
        mem=config["large_mem"],
    shell:
        " checkm lineage_wf "
        " --file {params.completeness} "
        " --tab_table "
        " --extension {params.extension} "
        " --threads {threads} "
        " {input.bins} "
        " {params.tmp_output_dir} &> {log} "
        " ;\n"
        " checkm tree_qa "
        " {params.tmp_output_dir} "
        " --out_format 2 "
        " --file {output.taxonomy} "
        " --tab_table "
        " &>> {log}"


## Combine bin stats
localrules:
    build_bin_report,
    combine_bin_stats,


if config["bin_quality_asesser"].lower() == "checkm2":

    rule combine_bin_stats:
        input:
            completeness_files=expand(
                "{sample}/binning/{{binner}}/bin_quality/checkm2/quality_report.tsv", sample=SAMPLES,
            ),
        output:
            bin_table="reports/genomic_bins_{binner}.tsv",
        params:
            samples=SAMPLES,
        log:
            "logs/binning/combine_stats_{binner}.log",
        script:
            "../scripts/combine_checkm2.py"


elif config["bin_quality_asesser"].lower() == "busco":

    rule combine_bin_stats:
        input:
            completeness_files=expand(
                "{sample}/binning/{{binner}}/bin_quality/busco.tsv", sample=SAMPLES
            ),
        output:
            bin_table="reports/genomic_bins_{binner}.tsv",
        params:
            samples=SAMPLES,
        log:
            "logs/binning/combine_stats_{binner}.log",
        script:
            "../scripts/combine_busco.py"


elif config["bin_quality_asesser"].lower() == "checkm":

    rule combine_bin_stats:
        input:
            completeness_files=expand(
                "{sample}/binning/{{binner}}/bin_quality/checkm.tsv", sample=SAMPLES
            ),
            taxonomy_files=expand(
                "{sample}/binning/{{binner}}/bin_quality/checkm_taxonomy.tsv",
                sample=SAMPLES,
            ),
        output:
            bin_table="reports/genomic_bins_{binner}.tsv",
        params:
            samples=SAMPLES,
        log:
            "logs/binning/combine_stats_{binner}.log",
        script:
            "../scripts/combine_checkm.py"


else:

    raise Exception(
        "'bin_quality_asesser' should be 'checkm2','busco' or 'checkm' got {bin_quality_asesser}".format(
            **config
        )
    )


rule build_bin_report:
    input:
        bin_table="reports/genomic_bins_{binner}.tsv",
    output:
        report="reports/bin_report_{binner}.html",
    conda:
        "%s/report.yaml" % CONDAENV
    log:
        "logs/binning/report_{binner}.log",
    script:
        "../report/bin_report.py"
