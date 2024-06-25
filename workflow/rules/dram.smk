DBDIR = config["database_dir"]


def get_dram_config(wildcards):
    old_dram_path = f"{DBDIR}/Dram"
    if Path(old_dram_path).exists():
        logger.error(
            f"Detected an old database for DRAM in {old_dram_path}. You can delete it."
        )

    return config.get("dram_config_file", f"{DBDIR}/DRAM/DRAM.config")


localrules:
    dram_download,
    concat_annotations,


rule dram_download:
    output:
        dbdir=directory(f"{DBDIR}/DRAM/db/"),
        config=f"{DBDIR}/DRAM/DRAM.config",
    threads: config["threads"]
    resources:
        mem_mb=config["mem"] * 1000,
        time_min=60 * config["runtime"]["default"],
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


rule DRAM_annotate:
    input:
        fasta="genomes/genomes/{genome}.fasta",
        #checkm= "genomes/checkm/completeness.tsv",
        #gtdb_dir= "genomes/taxonomy/gtdb/classify",
        config=get_dram_config,
    output:
        outdir=directory("genomes/annotations/dram/intermediate_files/{genome}"),
    threads: config["simplejob_threads"]
    resources:
        mem_mb=config["simplejob_mem"] * 1000,
        time_min=60 * config["runtime"]["default"],
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


rule concat_annotations:
    input:
        get_all_dram,
    output:
        expand("genomes/annotations/dram/{annotation}", annotation=DRAM_ANNOTATON_FILES),
    resources:
        time_min=60 * config["runtime"]["default"],
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
        config=get_dram_config,
    output:
        outdir=directory("genomes/annotations/dram/distil"),
    threads: 1
    resources:
        mem_mb=config["simplejob_mem"] * 1000,
        ttime_min=60 * config["runtime"]["simplejob"],
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
        annotations="genomes/annotations/dram/annotations.tsv",
        config=get_dram_config,
    output:
        "genomes/annotations/dram/kegg_modules.tsv",
    threads: 1
    resources:
        mem_mb=config["simplejob_mem"] * 1000,
        time_min=60 * config["runtime"]["default"],
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
