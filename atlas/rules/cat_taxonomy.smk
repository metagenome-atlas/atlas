


# this is a HACK because
localrules: get_genome_for_cat
rule get_genome_for_cat:
    input:
        "genomes/genomes"
    output:
        dynamic(temp("genomes/taxonomy/intermediate_files/{genome}/{genome}.fasta"))
    shadow:
        "shallow"
    run:

        import os,shutil
        genome_path= shutil.path.join(input[0],'{genome}.fasta')
        Genomes = glob_wildcards(genome_path).genome

        for genome in Genomes:
            os.makedirs(f"genomes/taxonomy/intermediate_files/{genome}")
            os.copy(genome_path.format(genome=genome), f"genomes/taxonomy/intermediate_files/{genome}/{genome}.fasta")



# CAT output files with 'CAT' as prefix
#CAT.bin2classification.txt  CAT.concatenated.alignment.diamond  CAT.concatenated.predicted_proteins.faa  CAT.log          summary.txt
#CAT.bin2name.txt            CAT.concatenated.fasta              CAT.concatenated.predicted_proteins.gff  CAT.ORF2LCA.txt

rule cat_on_bin:
    input:
        flag=CAT_flag_downloaded,
        genome= "genomes/taxonomy/intermediate_files/{genome}/{genome}.fasta",
        proteins= "genomes/annotations/genes/{genome}.faa"
    output:
        "genomes/taxonomy/intermediate_files/{genome}/{genome}.bin2classification.txt"
    params:
        db_folder=CAT_DIR,
        bin_folder=lambda wc,input: os.path.dirname(input.genome),
        extension=".fasta",
        out_prefix= lambda wc,output: os.path.join(os.path.dirname(output[0]),wc.genome)
    resources:
        mem= config['java_mem']
    threads:
        config['threads']
    conda:
        "%s/cat.yaml" % CONDAENV
    log:
        "logs/genomes/taxonomy/{genome}.log"
    shell:
        " CAT bins "
        " -b {params.bin_folder} -p {input.proteins} "
        "-d {params.db_folder} -t {params.db_folder} --nproc {threads} "
        " --bin_suffix {params.extension} "
        " --out_prefix {params.out_prefix} &> >(tee {log})"



localrules: merge_taxonomy, cat_get_name
rule merge_taxonomy:
    input:
        taxid=dynamic(rules.cat_on_bin.output),
    output:
        "genomes/taxonomy/taxonomy_ids.tsv"
    threads:
        1
    run:
        import pandas as pd
        out= pd.concat([pd.read_table(file,index_col=0) for file in input],axis=0).sort_index()

        out.to_csv(output[0],sep='\t')

rule cat_get_name:
    input:
        "genomes/taxonomy/taxonomy_ids.tsv"
    output:
        "genomes/taxonomy/taxonomy_names.tsv"
    params:
        db_folder=CAT_DIR,
    conda:
        "%s/cat.yaml" % CONDAENV
    threads:
        1
    shell:
        " CAT add_names -i {input} -t {params.db_folder} "
        " -o {output} --only_official "
