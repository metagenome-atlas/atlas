DBDIR = config["database_dir"]


def get_dram_config(wildcards):
    return config.get("dram_config_file", f"{DBDIR}/DRAM.config")


rule dram_download:
    output:
        dbdir=directory(f"{DBDIR}/Dram/"),
        config=f"{DBDIR}/DRAM.config",
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
        touch(f"{DBDIR}/dram_config_imported"),
    threads: 1
    conda:
        "../envs/dram.yaml"
    shell:
        "DRAM-setup.py import_config --config_loc {input}"


rule DRAM_annotate:
    input:
        fasta="genomes/genomes/{genome}.fasta",
        #checkm= "genomes/checkm/completeness.tsv",
        #gtdb_dir= "genomes/taxonomy/gtdb/classify",
        flag=rules.DRAM_set_db_loc.output,
    output:
        outdir=directory("genomes/annotations/dram/intermediate_files/{genome}"),
    threads: config["threads"]
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


DRAM_ANNOTATON_FILES = ["annotations.tsv", "rrnas.tsv", "trnas.tsv"]

DRAM_ANNOTATON_FILES = ["annotations.tsv", "rrnas.tsv", "trnas.tsv"]


localrules:
    concat_annotations,


rule concat_annotations:
    input:
        get_all_dram,
    output:
        expand("genomes/annotations/dram/{annotation}", annotation=DRAM_ANNOTATON_FILES),
    resources:
        time=config["runtime"]["default"],
    #     mem = config['mem']
    run:
        # from utils import io
        for i, annotation_file in enumerate(DRAM_ANNOTATON_FILES):

            input_files = [
                os.path.join(dram_folder, annotation_file) for dram_folder in input
            ]

            # drop files that don't exist for rrna and trna
            if not i == 0:
                input_files = [f for f in input_files if os.path.exists(f)]

            shell(f"head -n1 {input_files[0]} > {output[i]} ")
            for f in input_files:
                shell(f"tail -n+2 {f} >> {output[i]}")

        # io.pandas_concat(input_files, output[i],sep='\t',index_col=0, axis=0)



rule DRAM_destill:
    input:
        rules.concat_annotations.output,
        flag=rules.DRAM_set_db_loc.output,
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
        " --input_file {input[0]}"
        " --rrna_path {input[1]}"
        " --trna_path {input[2]}"
        " --output_dir {output} "
        "  &> {log}"


rule get_all_modules:
    input:
        "genomes/annotations/dram/annotations.tsv",
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
