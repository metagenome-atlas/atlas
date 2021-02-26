
DBDIR = config['database_dir']
genome_folder= os.path.dirname(config.get('genome_folder', 'genomes/genomes'))

def get_dram_config(wildcards):
    return config.get('dram_config_file', f"{DBDIR}/DRAM.config")


rule dram_download:
    output:
        dbdir= directory(f"{DBDIR}/Dram/"),
        config= f"{DBDIR}/DRAM.config"
    threads:
        config['threads']
    resources:
        mem= config['mem'],
        time= config["runtime"]["default"]
    log:
        "log/dram/download_dram.log"
    benchmark:
        "log/benchmarks/dram/download_dram.tsv"
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



localrules: DRAM_set_db_loc
rule DRAM_set_db_loc:
    input:
        get_dram_config
    output:
        touch("logs/dram_config_imported")
    threads:
        1
    conda:
        "../envs/dram.yaml"
    shell:
        "DRAM-setup.py import_config --config_loc {input}"

rule DRAM_annotate:
    input:
        fasta=f"{genome_folder}/{{genome}}.fasta",
        #checkm= "genomes/checkm/completeness.tsv",
        #gtdb_dir= "genomes/taxonomy/gtdb/classify",
        flag= rules.DRAM_set_db_loc.output
    output:
        "genomes/annotations/dram/intermediate_files/{genome}/annotations.tsv",
        "genomes/annotations/dram/intermediate_files/{genome}/trnas.tsv"
    threads:
        config['threads']
    resources:
        mem= config['simplejob_mem'],
        time= config['runtime']['default']
    conda:
        "../envs/dram.yaml"
    params:
        gtdb_file="gtdbtk.bac120.summary.tsv",
        outdir= "genomes/annotations/dram/intermediate_files/{genome}"
    log:
        "log/dram/run_dram/{genome}.log"
    benchmark:
        "log/benchmarks/dram/run_dram/{genome}.tsv"
    shell:
        "rm -r {params.outdir} ; "
        " DRAM.py annotate "
        " --input_fasta {input.fasta}"
        " --output_dir {params.outdir} "
        " --prodigal_mode single "
        #" --gtdb_taxonomy {input.gtdb_dir}/{params.gtdb_file} "
        #" --checkm_quality {input.checkm} "
        " --threads {threads} "
        " --verbose &> {log}"


def get_all_dram(wildcards):

    all_genomes = glob_wildcards(f"{genome_folder}/{{i}}.fasta").i

    return expand(rules.DRAM_annotate.output[0],
           genome=all_genomes)



rule dram:
    input:
        get_all_dram
