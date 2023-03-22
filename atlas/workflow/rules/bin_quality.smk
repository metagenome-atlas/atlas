bin_quality_input_folder = "{sample}/binning/{binner}/bins"


### Checkm2 ###




rule run_checkm2:
    input:
        fasta_dir=bin_quality_input_folder,
        db=rules.checkm2_download_db.output,
    output:
        dir=directory("{sample}/binning/{binner}/checkm2"),
        table = "{sample}/binning/{binner}/checkm2_report.tsv",
        faa=directory("{sample}/binning/{binner}/faa"),
    params:
        outdir=lambda wc, output: Path(output[0]).parent,
        lowmem=" --lowmem " if config["mem"] < 10 else "",
    conda:
        "../envs/checkm2.yaml"
    shadow:
        "minimal"
    threads: config["threads"]
    log:
        "{sample}/logs/binning/{binner}/checkm2.log",
        "{sample}/binning/{binner}/checkm2/checkm2.log",
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
        " --output-directory {output.dir} "
        " &> {log[0]} "
        ";\n"
        " cp {output.dir}/quality_report.tsv {output.table} 2>> {log[0]} ; "
        " mv {output.dir}/protein_files {output.faa} 2>> {log[0]} ; "




## GUNC ###


rule run_gunc:
    input:
        db=rules.download_gunc.output[0].format(**config),
        fasta_dir= "{sample}/binning/{binner}/faa"
    output:
        table="{sample}/binning/{binner}/gunc_output.tsv",
        folder=directory("{sample}/binning/{binner}/gunc")
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
        " gunc run "
        " --threads {threads} "
        " --gene_calls "
        " --db_file {input.db} "
        " --input_dir {input.fasta_dir} "
        " --temp_dir {resources.tmpdir} "
        " --file_suffix .fasta "
        " --out_dir {output.folder} &> {log} "
        " ;\n "
        " cp {output.folder}/*.tsv {output.table} 2>> {log}"


##### BUSCO  #########
'''



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

'''
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





rule combine_gunc:
    input:
        expand(
            "{sample}/binning/{{binner}}/gunc_output.tsv",
            sample=SAMPLES,
        ),
    output:
        bin_table="Binning/{binner}/gunc_report.tsv",
    params:
        samples=SAMPLES,
    log:
        "logs/binning/{binner}/combine_gunc.log",
    script:
        "../scripts/combine_gunc.py"



rule combine_bin_stats:
    input:
        completeness_files=expand(
            "{sample}/binning/{{binner}}/checkm2_report.tsv",
            sample=SAMPLES,
        )
    output:
        bin_table="Binning/{binner}/checkm2_report.tsv",
    params:
        samples=SAMPLES,
    log:
        "logs/binning/combine_stats_{binner}.log",
    script:
        "../scripts/combine_checkm2.py"


localrules:
    get_bin_filenames,

rule get_bin_filenames:
    input:
        dirs=expand(
            "{sample}/binning/{{binner}}/bins",
            sample=SAMPLES,
        ),
        protein_dirs=expand(
            "{sample}/binning/{{binner}}/faa",
            sample=SAMPLES,
        ),
    output:
        filenames="Binning/{binner}/paths.tsv",
    run:
        import pandas as pd
        from pathlib import Path
        from utils import io


        def get_list_of_files(dirs,pattern):

            fasta_files = []

            # searh for fasta files (.f*) in all bin folders
            for dir in dirs:
                dir = Path(dir)
                fasta_files += list(dir.glob(pattern))

            filenames = pd.DataFrame(fasta_files, columns=["Filename"])
            filenames.index = filenames.Filename.apply(io.simplify_path)
            filenames.index.name = "Bin"

            filenames.sort_index(inplace=True)

            return(filenames)

        fasta_filenames= get_list_of_files(input.dirs,"*.f*")
        faa_filenames= get_list_of_files(input.protein_dirs,"*.faa")

        assert faa_filenames.index == fasta_filenames.index, "faa index and faa index are nt the same"

        faa_filenames.name= "Proteins"

        filenames= pd.concat((fasta_filenames,faa_filenames),axis=1)

        filenames.to_csv(output.filenames, sep="\t")


rule calculate_stats:
    input:
        rules.get_bin_filenames.output.filenames,
    output:
        "Binning/{binner}/genome_stats.tsv",
    threads: config["threads"]
    run:
        from utils.genome_stats import get_many_genome_stats
        import pandas as pd

        filenames = pd.read_csv(input[0], sep="\t", index_col=0).squeeze()
        get_many_genome_stats(filenames, output[0], threads)

rule merge_bin_stats:
    input:
        stats= "Binning/{binner}/genome_stats.tsv",
        gunc= "Binning/{binner}/gunc_report.tsv",
        quality= "Binning/{binner}/checkm2_report.tsv"
    output:
        "Binning/{binner}/quality_report.tsv"



rule build_bin_report:
    input:
        bin_table="Binning/{binner}/quality_report.tsv",
    output:
        report="Binning/{binner}/report.html",
    conda:
        "%s/report.yaml" % CONDAENV
    log:
        "logs/binning/report_{binner}.log",
    script:
        "../report/bin_report.py"


