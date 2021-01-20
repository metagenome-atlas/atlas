
def get_dram_config(wildcards):

    if not 'dram_config_file' in config:

        raise IOError("If you want to run DRAM functional annotation. "
                      "You need to download the databases first "
                      "and give me the path to the DRAM config file.\n"
                      "Download the DRAM databases as described in the DRAM docs. "
                      "Then run:\n"
                      " DRAM-setup.py export_config --output_file /path/to/DRAM.config \n"
                      "Then set the path in the atlas config.yaml file:\n\n"
                      "dram_config_file: /path/to/DRAM.config"
                      )
    else:
        return config['dram_config_file']

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
        fasta="genomes/genomes/{genome}.fasta",
        checkm= "genomes/checkm/completeness.tsv",
        gtdb_dir= "genomes/taxonomy/gtdb/classify",
        flag= rules.DRAM_set_db_loc.output
    output:
        "genomes/annotations/dram/intermediate_files/{genome}/annotations.tsv",
        "genomes/annotations/dram/intermediate_files/{genome}/trnas.tsv"
    threads:
        config['threads']
    resources:
        mem=config['large_mem']
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
        "DRAM.py annotate "
        " --input_fasta {input.fasta}"
        " --output_dir {params.outdir} "
        " --gtdb_taxonomy {input.gtdb_dir}/{params.gtdb_file} "
        " --checkm_quality {input.checkm} "
        " --threads {threads} "
        " --verbose &> {log}"


def get_all_dram(wildcards):
    all_genomes = glob_wildcards("genomes/genomes/{i}.fasta").i
    return expand(rules.DRAM_annotate.output[0],
           genome=all_genomes)

rule dram:
    input:
        get_all_dram
