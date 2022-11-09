## dRep
localrules:
    all_contigs2bins,


rule all_contigs2bins:
    input:
        expand(
            "{sample}/binning/{binner}/cluster_attribution.tsv",
            sample=SAMPLES,
            binner=config["final_binner"],
        ),
    output:
        temp("genomes/clustering/all_contigs2bins.tsv.gz"),
    run:
        from utils.io import cat_files

        cat_files(input, output[0], gzip=True)


localrules:
    merge_checkm,
    filter_bins,


localrules:
    get_bin_filenames,


rule get_bin_filenames:
    input:
        dirs=expand(
            "{sample}/binning/{binner}/bins",
            sample=SAMPLES,
            binner=config["final_binner"],
        ),
    output:
        filenames="genomes/filter/paths.tsv",
    run:
        import pandas as pd
        from pathlib import Path
        from utils import io

        fasta_files = []


        # searh for fasta files (.f*) in all bin folders
        for dir in input.dirs:
            dir = Path(dir)
            fasta_files += list(dir.glob("*.f*"))

        filenames = pd.DataFrame(fasta_files, columns=["Filename"])
        filenames.index = filenames.Filename.apply(io.simplify_path)
        filenames.index.name = "Bin"

        filenames.to_csv(output.filenames, sep="\t")


rule calculate_stats:
    input:
        rules.get_bin_filenames.output.filenames,
    output:
        "genomes/all_bins/genome_stats.tsv",
    threads: config["threads"]
    run:
        from utils.genome_stats import get_many_genome_stats
        import pandas as pd

        filenames = pd.read_csv(input[0], sep="\t", index_col=0).squeeze()
        get_many_genome_stats(filenames, output[0], threads)


rule merge_checkm:
    input:
        completeness=expand(
            "{sample}/binning/{binner}/checkm/completeness.tsv",
            sample=SAMPLES,
            binner=config["final_binner"],
        ),
        taxonomy=expand(
            "{sample}/binning/{binner}/checkm/taxonomy.tsv",
            sample=SAMPLES,
            binner=config["final_binner"],
        ),
        markers=expand(
            "{sample}/binning/{binner}/checkm/storage/tree/concatenated.fasta",
            sample=SAMPLES,
            binner=config["final_binner"],
        ),
    output:
        checkm="genomes/all_bins/checkm_all_bins.tsv",
        markers="genomes/all_bins/all_bins_markers.fasta",
    run:
        import pandas as pd
        import shutil

        D = []

        for i in range(len(SAMPLES)):
            df = pd.read_csv(input.completeness[i], index_col=0, sep="\t")
            df = df.join(pd.read_csv(input.taxonomy[i], index_col=0, sep="\t"))
            D.append(df)

        D = pd.concat(D, axis=0)
        D.to_csv(output.checkm, sep="\t")

        with open(output.markers, "wb") as fout:
            for fasta in input.markers:
                shutil.copyfileobj(open(fasta, "rb"), fout)


rule filter_bins:
    input:
        paths=rules.get_bin_filenames.output.filenames,
        quality="genomes/all_bins/checkm_all_bins.tsv",
        stats="genomes/all_bins/genome_stats.tsv",
    output:
        quality="genomes/all_bins/filtered_quality.tsv",
        paths=temp("genomes/all_bins/filtered_bins.txt"),
        quality_for_derep=temp("genomes/all_bins/filtered_quality.csv"),
    threads: 1
    log:
        "logs/genomes/filter_bins.log",
    params:
        filter_criteria=config["genome_filter_criteria"],
    script:
        "../scripts/filter_genomes.py"


rule dereplication:
    input:
        paths="genomes/all_bins/filtered_bins.txt",
        quality="genomes/all_bins/filtered_quality.csv",
    output:
        dir=temp(directory("genomes/dereplicated_genomes")),
        mapping_file=temp("genomes/clustering/allbins2genome_oldname.tsv"),
    threads: config["threads"]
    log:
        "logs/genomes/dereplication.log",
    conda:
        "../envs/galah.yaml"
    params:
        ANI=config["genome_dereplication"]["ANI"],
        quality_formula=config["genome_dereplication"]["quality_formula"],
        opt_parameters=config["genome_dereplication"]["opt_parameters"],
        min_overlap=config["genome_dereplication"]["overlap"],
    shell:
        " galah cluster "
        " --genome-fasta-list {input.paths}"
        " --genome-info {input.quality} "
        " --ani {params.ANI} "
        " --min-aligned-fraction {params.min_overlap} "
        " {params.opt_parameters} "
        " --quality-formula {params.quality_formula} "
        " --threads {threads} "
        " --output-representative-fasta-directory {output.dir} "
        " --output-cluster-definition {output.mapping_file} "
        " --precluster-method finch "
        " &> {log} "


localrules:
    rename_genomes,


checkpoint rename_genomes:
    input:
        genomes="genomes/dereplicated_genomes",
        mapping_file="genomes/clustering/allbins2genome_oldname.tsv",
        genome_quality=f"reports/genomic_bins_{config['final_binner']}.tsv",
    output:
        dir=directory("genomes/genomes"),
        mapfile_contigs="genomes/clustering/contig2genome.tsv",
        mapfile_old2mag="genomes/clustering/old2newID.tsv",
        mapfile_allbins2mag="genomes/clustering/allbins2genome.tsv",
        genome_quality="genomes/genome_quality.tsv",
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
        mem=config["simplejob_mem"],
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

localrules: concat_orfs, concat_genomes

rule concat_orfs:
    input:
        lambda wc: get_all_genes(wc, extension=".fna")
    output:
        temp("genomes/all_orfs.fasta"),
    shell:
        "cat {input} > {output}"


rule concat_genomes:
    input:
        genome_dir,
    output:
        "genomes/all_contigs.fasta",
    params:
        ext="fasta",
    shell:
        "cat {input}/*{params.ext} > {output}"


rule index_genomes:
    input:
        target="genomes/all_contigs.fasta",
    output:
        "ref/genomes.mmi",
    log:
        "logs/genomes/alignment/index.log",
    params:
        index_size="12G",
    threads: 3
    resources:
        mem=config["mem"],
    wrapper:
        "v1.19.0/bio/minimap2/index"


rule align_reads_to_genomes:
    input:
        target=rules.index_genomes.output,
        query=get_quality_controlled_reads,
    output:
        "genomes/alignments/{sample}.bam",
    log:
        "logs/genomes/alignment/{sample}_map.log",
    params:
        extra="-x sr",
    threads: config["threads"]
    resources:
        mem=config["mem"],
        mem_mb=config["mem"] * 1000,
    wrapper:
        "v1.19.0/bio/minimap2/aligner"



rule pileup_MAGs:
    input:
        bam="genomes/alignments/{sample}.bam",
        orf= "genomes/all_orfs.fasta"
    output:
        covstats=temp("genomes/alignments/coverage/{sample}.tsv.gz"),
        bincov=temp("genomes/alignments/coverage_binned/{sample}.tsv.gz"),
        orf= "genomes/alignments/orf_coverage/{sample}.tsv.gz"
    log:
        "logs/genomes/alignments/pilup_{sample}.log",
    conda:
        "../envs/required_packages.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
    shell:
        "pileup.sh in={input.bam} "
        " threads={threads} "
        " -Xmx{resources.java_mem}G "
        " covstats={output.covstats} "
        " fastaorf={input.orf} outorf={output.orf} "
        " concise=t "
        " physical=t "
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
        coverage_contigs = "genomes/counts/coverage_contigs.parquet",
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


# TODO mapping rate from pileup
