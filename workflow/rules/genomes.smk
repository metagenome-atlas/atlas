

localrules:
    rename_genomes,


checkpoint rename_genomes:
    input:
        paths="Binning/{binner}/raw_bins/paths.tsv".format(
            binner=config["final_binner"]
        ),
        mapping_file="Binning/{binner}/bins2species.tsv".format(
            binner=config["final_binner"]
        ),
        genome_info=f"Binning/{config['final_binner']}/filtered_bin_info.tsv",
    output:
        dir=directory("genomes/genomes"),
        mapfile_contigs="genomes/clustering/contig2genome.tsv",
        mapfile_old2mag="genomes/clustering/old2newID.tsv",
        mapfile_allbins2mag="genomes/clustering/allbins2genome.tsv",
        genome_info="genomes/genome_quality.tsv",
    params:
        rename_contigs=config["rename_mags_contigs"],
    shadow:
        "shallow"
    log:
        "logs/genomes/rename_genomes.log",
    script:
        "../scripts/rename_genomes.py"


def get_genome_dir():
    if ("genome_dir" in config) and (config["genome_dir"] is not None):
        genome_dir = config["genome_dir"]
        assert os.path.exists(genome_dir), f"{genome_dir} Doesn't exists"

        logger.info(f"Set genomes from {genome_dir}.")

        # check if genomes are present
        genomes = glob_wildcards(os.path.join(genome_dir, "{genome}.fasta")).genome

        if len(genomes) == 0:
            logger.error(f"No genomes found with fasta extension in {genome_dir} ")
            exit(1)

    else:
        genome_dir = "genomes/genomes"

    return genome_dir


genome_dir = get_genome_dir()


def get_all_genomes(wildcards):
    global genome_dir

    if genome_dir == "genomes/genomes":
        checkpoints.rename_genomes.get()

    # check if genomes are present
    genomes = glob_wildcards(os.path.join(genome_dir, "{genome}.fasta")).genome

    if len(genomes) == 0:
        logger.error(
            f"No genomes found with fasta extension in {genome_dir} "
            "You don't have any Metagenome assembled genomes with sufficient quality. "
            "You may want to change the assembly, binning or filtering parameters. "
            "Or focus on the genecatalog workflow only."
        )
        exit(1)

    return genomes


rule get_contig2genomes:
    input:
        genome_dir,
    output:
        "genomes/clustering/contig2genome.tsv",
    run:
        from glob import glob

        fasta_files = glob(input[0] + "/*.f*")

        with open(output[0], "w") as out_contigs:
            for fasta in fasta_files:
                bin_name, ext = os.path.splitext(os.path.split(fasta)[-1])
                # if gz remove also fasta extension
                if ext == ".gz":
                    bin_name = os.path.splitext(bin_name)[0]

                    # write names of contigs in mapping file
                with open(fasta) as f:
                    for line in f:
                        if line[0] == ">":
                            header = line[1:].strip().split()[0]
                            out_contigs.write(f"{header}\t{bin_name}\n")


# alternative way to get to contigs2genomes for quantification with external genomes


ruleorder: get_contig2genomes > rename_genomes


# rule predict_genes_genomes:
#     input:
#         dir= genomes_dir
#     output:
#         directory("genomes/annotations/genes")
#     conda:
#         "%s/prodigal.yaml" % CONDAENV
#     log:
#         "logs/genomes/prodigal.log"
#     shadow:
#         "shallow"
#     threads:
#         config.get("threads", 1)
#     script:
#         "predict_genes_of_genomes.py"


rule predict_genes_genomes:
    input:
        os.path.join(genome_dir, "{genome}.fasta"),
    output:
        fna="genomes/annotations/genes/{genome}.fna",
        faa="genomes/annotations/genes/{genome}.faa",
        gff=temp("genomes/annotations/genes/{genome}.gff"),
    conda:
        "%s/prodigal.yaml" % CONDAENV
    log:
        "logs/genomes/prodigal/{genome}.txt",
    threads: 1
    resources:
        mem_mb=config["simplejob_mem"] * 1024,
        time=config["runtime"]["simplejob"],
    shell:
        """
        prodigal -i {input} -o {output.gff} -d {output.fna} \
            -a {output.faa} -p meta -f gff 2> {log}
        """


def get_all_genes(wildcards, extension=".faa"):
    return expand(
        "genomes/annotations/genes/{genome}{extension}",
        genome=get_all_genomes(wildcards),
        extension=extension,
    )


localrules:
    all_prodigal,


rule all_prodigal:
    input:
        get_all_genes,
    output:
        touch("genomes/annotations/genes/predicted"),


### Quantification


localrules:
    concat_orfs,
    concat_genomes,


rule concat_orfs:
    input:
        lambda wc: get_all_genes(wc, extension=".fna"),
    output:
        temp("genomes/all_orfs.fasta"),
    shell:
        "cat {input} > {output}"


rule concat_genomes:
    input:
        genome_dir,
    output:
        temp("genomes/all_contigs.fasta"),
    params:
        ext="fasta",
    shell:
        "cat {input}/*{params.ext} > {output}"


if config["genome_aligner"] == "minimap":

    rule index_genomes:
        input:
            target=ancient("genomes/all_contigs.fasta"),
            timestamp=genome_dir,
        output:
            "ref/genomes.mmi",
        log:
            "logs/genomes/alignmentsindex.log",
        params:
            index_size="12G",
        threads: 3
        resources:
            mem_mb=config["mem"] * 1024,
        wrapper:
            "v1.19.0/bio/minimap2/index"

    rule align_reads_to_genomes:
        input:
            target=rules.index_genomes.output,
            query=get_quality_controlled_reads,
        output:
            "genomes/alignments/bams/{sample}.bam",
        log:
            "logs/genomes/alignments/{sample}_map.log",
        params:
            extra="-x sr",
            sort="coordinate",
        threads: config["threads"]
        resources:
            mem_mb=config["mem"] * 1024,
            mem_mb=config["mem"] * 1000,
        wrapper:
            "v1.19.0/bio/minimap2/aligner"

elif config["genome_aligner"] == "bwa":

    rule index_genomes:
        input:
            ancient("genomes/all_contigs.fasta"),
            timestamp=genome_dir,
        output:
            multiext("ref/genomes", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
        log:
            "logs/genomes/alignments/bwa_index.log",
        threads: 4
        resources:
            mem_mb=config["mem"] * 1024,
        wrapper:
            "v1.19.0/bio/bwa-mem2/index"

    rule align_reads_to_genomes:
        input:
            idx=rules.index_genomes.output,
            reads=get_quality_controlled_reads,
        output:
            "genomes/alignments/bams/{sample}.bam",
        log:
            "logs/genomes/alignments/{sample}_bwa.log",
        params:
            extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
            sort="samtools",
            sort_order="coordinate",
        threads: config["threads"]
        resources:
            mem_mb=config["mem"] * 1024,
            mem_mb=config["mem"] * 1000,
        wrapper:
            "v1.19.0/bio/bwa-mem2/mem"

else:
    raise Exception(
        "'genome_aligner' not understood, it should be 'minimap' or 'bwa', got '{genome_aligner}'. check config file".format(
            **config
        )
    )


# path change for bam file
localrules:
    move_old_bam,


ruleorder: move_old_bam > align_reads_to_genomes


rule move_old_bam:
    input:
        "genomes/alignments/{sample}.bam",
    output:
        "genomes/alignments/bams/{sample}.bam",
    log:
        "logs/genomes/alignments/{sample}_move.log",
    shell:
        "mv {input} {output} > {log}"


rule mapping_stats_genomes:
    input:
        bam="genomes/alignments/bams/{sample}.bam",
    output:
        "genomes/alignments/stats/{sample}.stats",
    log:
        "logs/genomes/alignments/{sample}_stats.log",
    threads: 1
    resources:
        mem_mb=config["simplejob_mem"] * 1024,
    wrapper:
        "v1.19.0/bio/samtools/stats"


rule multiqc_mapping_genome:
    input:
        expand("genomes/alignments/stats/{sample}.stats", sample=SAMPLES),
    output:
        "reports/genome_mapping/results.html",
    log:
        "logs/genomes/alignment/multiqc.log",
    wrapper:
        "v3.3.6/bio/multiqc"


rule pileup_MAGs:
    input:
        bam="genomes/alignments/bams/{sample}.bam",
        orf="genomes/all_orfs.fasta",
    output:
        covstats=temp("genomes/alignments/coverage/{sample}.tsv.gz"),
        bincov=temp("genomes/alignments/coverage_binned/{sample}.tsv.gz"),
        orf="genomes/alignments/orf_coverage/{sample}.tsv.gz",
    params:
        minmapq=config["minimum_map_quality"],
    log:
        "logs/genomes/alignments/pilup_{sample}.log",
    conda:
        "../envs/required_packages.yaml"
    threads: config["threads"]
    resources:
        mem_mb=config["mem"] * 1024,
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
    shell:
        "pileup.sh in={input.bam} "
        " threads={threads} "
        " -Xmx{resources.java_mem}G "
        " covstats={output.covstats} "
        " fastaorf={input.orf} outorf={output.orf} "
        " concise=t "
        " physical=t "
        " minmapq={params.minmapq} "
        " bincov={output.bincov} "
        " 2> {log}"


rule combine_coverages_MAGs:
    input:
        binned_coverage_files=expand(
            "genomes/alignments/coverage_binned/{sample}.tsv.gz", sample=SAMPLES
        ),
        coverage_files=expand(
            "genomes/alignments/coverage/{sample}.tsv.gz", sample=SAMPLES
        ),
        contig2genome="genomes/clustering/contig2genome.tsv",
    params:
        samples=SAMPLES,
    output:
        coverage_contigs="genomes/counts/coverage_contigs.parquet",
        counts="genomes/counts/counts_genomes.parquet",
        binned_cov="genomes/counts/binned_coverage.parquet",
        median_abund="genomes/counts/median_coverage_genomes.parquet",
    log:
        "logs/genomes/counts/combine_binned_coverages_MAGs.log",
    threads: 1
    resources:
        mem_mb=1000 * config["simplejob_mem"],
        time_min=config["runtime"]["simplejob"] * 60,
    script:
        "../scripts/combine_coverage_MAGs.py"
