# copy genome so CAt has only one genome per bin folder and paralelizes
localrules:
    get_genome_for_cat,


rule get_genome_for_cat:
    input:
        genomes="genomes/genomes/{genome}.fasta",
    output:
        genomes=temp("genomes/taxonomy/{genome}/{genome}.fasta"),
    run:
        import os, shutil

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        shutil.copy(input[0], output[0])


# CAT output files with 'CAT' as prefix
# CAT.bin2classification.txt  CAT.concatenated.alignment.diamond  CAT.concatenated.predicted_proteins.faa  CAT.log          summary.txt
# CAT.bin2name.txt            CAT.concatenated.fasta              CAT.concatenated.predicted_proteins.gff  CAT.ORF2LCA.txt


rule cat_on_bin:
    input:
        flag=CAT_flag_downloaded,
        genome="genomes/taxonomy/{genome}/{genome}.fasta",
    output:
        expand(
            "genomes/taxonomy/{{genome}}/{{genome}}.{extension}",
            extension=[
                "bin2classification.txt",
                "concatenated.alignment.diamond",
                "ORF2LCA.txt",
                "log",
            ],
        ),
    params:
        db_folder=CAT_DIR,
        bin_folder=lambda wc, input: os.path.dirname(input.genome),
        extension=".fasta",
        out_prefix=lambda wc, output: os.path.join(
            os.path.dirname(output[0]), wc.genome
        ),
        r=config["cat_range"],
        f=config["cat_fraction"],
    resources:
        mem=config["mem"],
    threads: config["threads"]
    conda:
        "%s/cat.yaml" % CONDAENV
    shadow:
        "shallow"
    log:
        "logs/genomes/taxonomy/{genome}.log",
    shell:
        " CAT bins "
        " -b {params.bin_folder} "
        " -f {params.f} -r {params.r} "
        " -d {params.db_folder} -t {params.db_folder} --nproc {threads} "
        " --bin_suffix {params.extension} "
        " --out_prefix {params.out_prefix} &> {log}"


def merge_taxonomy_input(wildcards):
    genomes = get_genomes_(wildcards)
    return expand(
        "genomes/taxonomy/{genome}/{genome}.bin2classification.txt", genome=genomes
    )


localrules:
    merge_taxonomy,
    cat_get_name,


rule merge_taxonomy:
    input:
        merge_taxonomy_input,
    output:
        "genomes/taxonomy/taxonomy_ids.tsv",
    threads: 1
    run:
        import pandas as pd

        out = pd.concat(
            [pd.read_csv(file, index_col=0, sep="\t") for file in input], axis=0
        ).sort_index()

        out.to_csv(output[0], sep="\t")


localrules:
    cat_get_name,
    parse_cat_output,


rule cat_get_name:
    input:
        "genomes/taxonomy/taxonomy_ids.tsv",
    output:
        "genomes/taxonomy/taxonomy_names.tsv",
    params:
        db_folder=CAT_DIR,
    conda:
        "%s/cat.yaml" % CONDAENV
    threads: 1
    shell:
        " CAT add_names -i {input} -t {params.db_folder} "
        " -o {output} --only_official "


rule parse_cat_output:
    input:
        rules.cat_get_name.output[0],
    output:
        "genomes/taxonomy/taxonomy.tsv",
    threads: 1
    script:
        "../scripts/parse_cat_taxonomy.py"
