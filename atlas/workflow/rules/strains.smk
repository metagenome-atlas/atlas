
rule instrain_profile:
    input:
        sam="genomes/alignments/{sample}.sam",
        genomes="genomes/all_contigs.fasta",
        genes=lambda wc: get_all_genes(wc, extension=".fna"),
        scaffold_to_genome="genomes/clustering/contig2genome.tsv",
    output:
        directory("strains/intermediate_files/{sample}"),
    threads: config["threads"]
    params:
        extra=config.get("instrain_profile_extra", ""),
    log:
        "logs/strains/profile/{sample}.log",
    conda:
        "../envs/instrain.yaml"
    benchmark:
        "logs/benchmarks/strains/profile/{sample}.tsv"
    shell:
        " cat {input.genes} > {resources.tmpdir}/all_genome_genes.fna 2> {log} "
        " ; "
        "inStrain profile "
        " {input.sam} {input.genomes} "
        " -o {output} "
        " -p {threads} "
        " -g {resources.tmpdir}/all_genome_genes.fna "
        " -s {input.scaffold_to_genome} "
        " --database_mode "
        " {params.extra} &>> {log}"


rule instrain_compare:
    input:
        profiles=expand("strains/intermediate_files/{sample}", sample=SAMPLES),
        scaffold_to_genome="genomes/clustering/contig2genome.tsv",
    output:
        directory("strains/comparison"),
    threads: config["threads"]
    params:
        extra=config.get("instrain_compare_extra", ""),
    log:
        "logs/strains/compare.log",
    conda:
        "../envs/instrain.yaml"
    benchmark:
        "logs/benchmarks/strains/compare.tsv"
    shell:
        "inStrain compare "
        " --input {input.profiles} "
        " -o {output} "
        " -p {threads} "
        " -s {input.scaffold_to_genome} "
        " --database_mode "
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


rule strains:
    input:
        "strains/comparison",
