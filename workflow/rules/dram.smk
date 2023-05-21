DBDIR = config["database_dir"]


def get_dram_config(wildcards):
    old_dram_path = f"{DBDIR}/Dram"
    if Path(old_dram_path).exists():
        logger.error(
            f"Detected an old database for DRAM in {old_dram_path}. You can delete it."
        )

    return config.get("dram_config_file", f"{DBDIR}/DRAM/DRAM.config")


rule dram_download:
    output:
        dbdir=directory(f"{DBDIR}/DRAM/db/"),
        config=f"{DBDIR}/DRAM/DRAM.config",
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["default"],
    log:
        "logs/dram/download_dram.log",
    benchmark:
        "logs/benchmarks/dram/download_dram.tsv"
    conda:
        "../envs/dram.yaml"
    shell:
        " DRAM-setup.py prepare_databases "
        " --output_dir {output.dbdir} "
        " --threads {threads} "
        " --verbose "
        " --skip_uniref "
        " &> {log} "
        " ; "
        " DRAM-setup.py export_config --output_file {output.config}"


localrules:
    DRAM_set_db_loc,


rule DRAM_set_db_loc:
    input:
        get_dram_config,
    output:
        touch(f"{DBDIR}/DRAM/dram_config_imported"),
    threads: 1
    conda:
        "../envs/dram.yaml"
    log:
        "logs/dram/set_db_loc.log",
    shell:
        "DRAM-setup.py import_config --config_loc {input} &> {log}"


rule DRAM_annotate:
    input:
        fasta="genomes/genomes/{genome}.fasta",
        #checkm= "genomes/checkm/completeness.tsv",
        #gtdb_dir= "genomes/taxonomy/gtdb/classify",
        flag=rules.DRAM_set_db_loc.output,
        config=f"{DBDIR}/DRAM/DRAM.config",
    output:
        outdir=directory("genomes/annotations/dram/intermediate_files/{genome}"),
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["default"],
    conda:
        "../envs/dram.yaml"
    params:
        extra=config.get("dram_extra", ""),
        min_contig_size=config.get("minimum_contig_length", "1000"),
    log:
        "logs/dram/run_dram/{genome}.log",
    benchmark:
        "logs/benchmarks/dram/run_dram/{genome}.tsv"
    shell:
        " DRAM.py annotate "
        " --config_loc {input.config} "
        " --input_fasta {input.fasta}"
        " --output_dir {output.outdir} "
        " --threads {threads} "
        " --min_contig_size {params.min_contig_size} "
        " {params.extra} "
        " --verbose &> {log}"
        #" --gtdb_taxonomy {input.gtdb_dir}/{params.gtdb_file} "
        #" --checkm_quality {input.checkm} "


def get_all_dram(wildcards):
    all_genomes = get_all_genomes(wildcards)

    return expand(rules.DRAM_annotate.output.outdir, genome=all_genomes)


DRAM_ANNOTATON_FILES = ["annotations.tsv"]


localrules:
    concat_annotations,


rule concat_annotations:
    input:
        get_all_dram,
    output:
        expand("genomes/annotations/dram/{annotation}", annotation=DRAM_ANNOTATON_FILES),
    resources:
        time=config["runtime"]["default"],
    run:
        from utils import io

        for i, annotation_file in enumerate(DRAM_ANNOTATON_FILES):
            input_files = [
                os.path.join(dram_folder, annotation_file) for dram_folder in input
            ]

            io.pandas_concat(
                input_files, output[i], sep="\t", index_col=0, axis=0, disk_based=True
            )


rule DRAM_destill:
    input:
        rules.concat_annotations.output,
        flag=rules.DRAM_set_db_loc.output,
        config=f"{DBDIR}/DRAM/DRAM.config",
    output:
        outdir=directory("genomes/annotations/dram/distil"),
    threads: 1
    resources:
        mem=config["simplejob_mem"],
        ttime=config["runtime"]["simplejob"],
    conda:
        "../envs/dram.yaml"
    log:
        "logs/dram/distil.log",
    shell:
        " DRAM.py distill "
        " --config_loc {input.config} "
        " --input_file {input[0]}"
        " --output_dir {output} "
        "  &> {log}"


rule get_all_modules:
    input:
        "genomes/annotations/dram/annotations.tsv",
        f"{DBDIR}/DRAM/DRAM.config",
    output:
        "genomes/annotations/dram/kegg_modules.tsv",
    threads: 1
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["default"],
    conda:
        "../envs/dram.yaml"
    log:
        "logs/dram/get_all_modules.log",
    script:
        "../scripts/DRAM_get_all_modules.py"


rule dram:
    input:
        "genomes/annotations/dram/distil",
        "genomes/annotations/dram/kegg_modules.tsv",
