bin_quality_input_folder = "{sample}/binning/{binner}/bins"




rule calculate_stats:
    input:
        bin_quality_input_folder,
    output:
        "{sample}/binning/{binner}/genome_stats.tsv",
    threads: config["threads"]
    params:
        extension=".fasta"
    run:
        from utils.genome_stats import get_many_genome_stats

        filenames = list(Path(input[0]).glob("*"+params.extension))

        get_many_genome_stats(filenames, output[0], threads)



rule combine_bin_stats:
    input:
        expand(
            "{sample}/binning/{{binner}}/genome_stats.tsv",
            sample=SAMPLES,
        ),
    output:
        "Binning/{binner}/genome_stats.tsv"
    params:
        samples=SAMPLES,
    log:
        "logs/binning/{binner}/combine_stats.log",
    run:

        try:
            from utils.io import pandas_concat

            pandas_concat(
            input,
            output[0]
            )

        except Exception as e:

            import traceback
            with open(log[0],"w") as logfile:
                traceback.print_exc(file=logfile)

            raise e
















### Checkm2 ###




rule run_checkm2:
    input:
        fasta_dir=bin_quality_input_folder,
        db=rules.checkm2_download_db.output,
    output:
        table = "{sample}/binning/{binner}/checkm2_report.tsv",
        faa=directory("{sample}/binning/{binner}/faa"),
    params:
        lowmem=" --lowmem " if config["mem"] < 10 else "",
        dir = lambda wc, output: Path(output.table).parent/"checkm2"
    conda:
        "../envs/checkm2.yaml"
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
        " --output-directory {params.dir} "
        " &> {log[0]} "
        ";\n"
        " cp {params.dir}/quality_report.tsv {output.table} 2>> {log[0]} ; "
        " mv {params.dir}/protein_files {output.faa} 2>> {log[0]} ; "





## GUNC ###


rule run_gunc:
    input:
        db=rules.download_gunc.output[0].format(**config),
        fasta_dir= "{sample}/binning/{binner}/faa"
    output:
        table="{sample}/binning/{binner}/gunc_output.tsv",
        folder=directory("{sample}/binning/{binner}/gunc"),
    params:
        extension=".faa"
    conda:
        "../envs/gunc.yaml"
    threads: config["threads"]
    log:
        "{sample}/logs/binning/{binner}/gunc.log",
    resources:
        time=int(config["runtime"]["default"]),
        mem_mb=config["mem"],
    shell:
        " mkdir {output.folder} 2> {log}"
        " ;\n"
        " gunc run "
        " --threads {threads} "
        " --gene_calls "
        " --db_file {input.db} "
        " --input_dir {input.fasta_dir} "
        " --temp_dir {resources.tmpdir} "
        " --file_suffix {params.extension} "
        " --out_dir {output.folder} &>> {log} "
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


## Combine bin stats
localrules:
    build_bin_report,
    combine_checkm2,
    combine_gunc





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
    run:

        try:
            from utils.io import pandas_concat

            pandas_concat(
            input,
            output[0]
            )

        except Exception as e:

            import traceback
            with open(log[0],"w") as logfile:
                traceback.print_exc(file=logfile)

            raise e




rule combine_checkm2:
    input:
        completeness_files=expand(
            "{sample}/binning/{{binner}}/checkm2_report.tsv",
            sample=SAMPLES,
        )
    output:
        bin_table="Binning/{binner}/checkm2_quality_report.tsv",
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

        assert all(faa_filenames.index == fasta_filenames.index), "faa index and faa index are nt the same"

        faa_filenames.columns= ["Proteins"]

        filenames= pd.concat((fasta_filenames,faa_filenames),axis=1)

        filenames.to_csv(output.filenames, sep="\t")

'''
rule merge_bin_info:
    input:
        stats ="Binning/{binner}/genome_stats.tsv",
        gunc= "Binning/{binner}/gunc_report.tsv",
        quality= "Binning/{binner}/checkm2_quality_report.tsv"
    output:
        "Binning/{binner}/combined_bin_info.tsv"

'''

rule build_bin_report:
    input:
        bin_table = "Binning/{binner}/checkm2_quality_report.tsv"
    output:
        report="Binning/{binner}/report.html",
    conda:
        "../envs/report.yaml"
    log:
        "logs/binning/report_{binner}.log",
    script:
        "../report/bin_report.py"



localrules:
    all_contigs2bins,


rule all_contigs2bins:
    input:
        expand(
            "{sample}/binning/{{binner}}/cluster_attribution.tsv",
            sample=SAMPLES,
        ),
    output:
        temp("Binning/{binner}/contigs2bins.tsv.gz"),
    run:
        from utils.io import cat_files

        cat_files(input, output[0], gzip=True)





rule quality_filter_bins:
    input:
        paths=rules.get_bin_filenames.output.filenames,
        stats ="Binning/{binner}/genome_stats.tsv",
        gunc= "Binning/{binner}/gunc_report.tsv",
        quality= "Binning/{binner}/checkm2_quality_report.tsv"
    output:
        info="Binning/{binner}/filtered_bin_info.tsv",
        paths=temp("Binning/{binner}/filtered_bins_paths.txt"),
        quality_for_derep=temp("Binning/{binner}/filtered_quality.csv"),
    threads: 1
    log:
        "logs/Binning/{binner}/filter_bins.log",
    params:
        filter_criteria=config["genome_filter_criteria"],
    script:
        "../scripts/filter_genomes.py"
