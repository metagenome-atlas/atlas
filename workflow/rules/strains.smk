
localrules:
    concatenate_genes_for_instrain,


rule concatenate_genes_for_instrain:
    input:
        lambda wc: get_all_genes(wc, extension=".fna"),
    output:
        "Intermediate/all_genome_genes.fna",
    log:
        "logs/strains/concatenate_genes.log",
    shell:
        "cat {input} > {output} 2> {log}"


rule instrain_profile:
    input:
        bam="genomes/alignments/bams/{sample}.bam",
        genomes="genomes/all_contigs.fasta",
        genes="Intermediate/all_genome_genes.fna",
        scaffold_to_genome="genomes/clustering/contig2genome.tsv",
    output:
        directory("Intermediate/strains/profiles/{sample}"),
    threads: config["threads"]
    params:
        extra=config.get("instrain_profile_extra", ""),
    log:
        "logs/strains/profile/{sample}.log",
    conda:
        "../envs/instrain.yaml"
    benchmark:
        "logs/benchmarks/strains/profile/{sample}.tsv"
    resources:
        mem_mb=config["large_mem"] * 1000,
        time_min=60 * config["runtime"]["long"],
    shell:
        "inStrain profile "
        " {input.bam} {input.genomes} "
        " -o {output} "
        " -p {threads} "
        " -s {input.scaffold_to_genome} "
        " --database_mode "
        " -g {input.genes} "
        " {params.extra} &>> {log}"


rule instrain_compare_genome:
    input:
        profiles=expand("Intermediate/strains/profiles/{sample}", sample=SAMPLES),
        scaffold_to_genome="genomes/clustering/contig2genome.tsv",
        genome=f"{genome_dir}/{{genome}}.fasta",
    output:
        directory("strains/comparison/{genome}"),
    threads: config["threads"]
    params:
        extra=config.get("instrain_compare_extra", ""),
    log:
        "logs/strains/compare/{genome}.log",
    conda:
        "../envs/instrain.yaml"
    benchmark:
        "logs/benchmarks/strains/compare_{genome}.tsv"
    resources:
        mem_mb=config["large_mem"] * 1000,
        time_min=60 * config["runtime"]["default"],
    shell:
        "inStrain compare "
        " --input {input.profiles} "
        " --stb {input.genome} "
        " -o {output} "
        " -p {threads} "
        " {params.extra} &> {log}"


# usage: inStrain compare -i [INPUT [INPUT ...]] [-o OUTPUT] [-p PROCESSES] [-d]
#                         [-h] [--version] [-s [STB [STB ...]]] [-c MIN_COV]
#                         [-f MIN_FREQ] [-fdr FDR] [--database_mode]
#                         [--breadth BREADTH] [-sc SCAFFOLDS] [--genome GENOME]
#                         [--store_coverage_overlap]
#                         [--store_mismatch_locations]
#                         [--include_self_comparisons] [--skip_plot_generation]
#                         [--group_length GROUP_LENGTH] [--force_compress]
#                         [-ani ANI_THRESHOLD] [-cov COVERAGE_TRESHOLD]
#                         [--clusterAlg {ward,single,complete,average,weighted,median,centroid}]


rule all_instrain_compare:
    input:
        lambda wildcards: expand(
            "strains/comparison/{genome}", genome=get_all_genomes(wildcards)
        ),
