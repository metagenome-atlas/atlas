

rule concat_genomes:
    input:
        genome_dir,
    output:
        "genomes/all_contigs.fasta"
    params:
        ext="fasta"
    shell:
        "cat {input}/*{params.ext} > {output}"

rule index_genomes:
    input:
        "genomes/all_contigs.fasta",
    output:
        "ref/genomes.mmi",
    log:
        "logs/genomes/alignment/index.log",
    params:
        index_size="12G",
    threads: 3
    resources:
        mem=config["mem"],
    conda:
        "../envs/minimap.yaml"
    shell:
        " minimap2 "
        " -I {params.index_size} "
        " -t {threads} "
        " -d {output} "
        " {input} 2> {log} "


rule align_reads_to_genomes:
    input:
        mmi=rules.index_genomes.output,
        reads=get_quality_controlled_reads,
    output:
        "genomes/alignments/{sample}.bam",
    log:
        "logs/genomes/alignment/{sample}_map.log",
    threads: config["threads"]
    resources:
        mem=config["mem"],
        mem_mb=config["mem"] * 1000,
    conda:
        "../envs/minimap.yaml"
    shell:
        " ( minimap2 "
        " -t {threads} "
        "-a "
        "-x sr "
        " {input.mmi} "
        " - "
        " | samtools view -b "
        " > {output} ) 2>{log}"



rule instrain_profile:
    input:
        sam="genomes/alignments/{sample}.bam",
        genomes="genomes/all_contigs.fasta",
        # genes=lambda wc: get_all_genes(wc, extension=".fna"),
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
    resources:
        mem=config["mem"],
        time=config["runtime"]["long"],
    shell:
        #" cat {input.genes} > {resources.tmpdir}/all_genome_genes.fna 2> {log} "
        #" ; "
        "inStrain profile "
        " {input.sam} {input.genomes} "
        " -o {output} "
        " -p {threads} "

        " -s {input.scaffold_to_genome} "
        " --database_mode "
        " {params.extra} &>> {log}"
        #" -g {resources.tmpdir}/all_genome_genes.fna "


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
    resources:
        mem=config["mem"],
        time=config["runtime"]["long"],
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
