## dRep
localrules:
    get_all_bins,
    all_contigs2bins,


rule get_all_bins:
    input:
        bins=expand(
            "{sample}/binning/{binner}/bins",
            sample=SAMPLES,
            binner=config["final_binner"],
        ),
    output:
        temp(directory("genomes/all_bins")),
    run:
        os.mkdir(output[0])
        from utils.io import symlink_relative

        for bin_folder in input.bins:

            fasta_files = [f for f in os.listdir(bin_folder) if f.endswith(".fasta")]
            symlink_relative(fasta_files, bin_folder, output[0])


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
    get_quality_for_dRep_from_checkm,
    merge_checkm,


rule get_quality_for_dRep_from_checkm:
    input:
        "reports/genomic_bins_{binner}.tsv".format(binner=config["final_binner"]),
    output:
        temp("genomes/quality.csv"),
    run:
        import pandas as pd

        D = pd.read_csv(input[0], index_col=0, sep="\t")

        D.index = D.index.astype(str) + ".fasta"
        D.index.name = "genome"
        D.columns = D.columns.str.lower()
        D.iloc[:, :3].to_csv(output[0])


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
        checkm="genomes/checkm/checkm_all_bins.tsv",
        markers="genomes/checkm/all_bins_markers.fasta",
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


rule dereplication:
    input:
        "genomes/all_bins",
        quality="genomes/quality.csv",
    output:
        directory("genomes/Dereplication/dereplicated_genomes"),
    threads: config["threads"]
    log:
        "logs/genomes/dereplication.log",
    shadow:
        "shallow"
    conda:
        "%s/dRep.yaml" % CONDAENV
    params:
        filter=(
            " --noQualityFiltering "
            if config["genome_dereplication"]["filter"]["noFilter"]
            else ""
        ),
        filter_length=config["genome_dereplication"]["filter"]["length"],
        filter_completeness=config["genome_dereplication"]["filter"]["completeness"],
        filter_contamination=config["genome_dereplication"]["filter"]["contamination"],
        ANI=config["genome_dereplication"]["ANI"],
        overlap=config["genome_dereplication"]["overlap"],
        completeness_weight=config["genome_dereplication"]["score"]["completeness"],
        contamination_weight=config["genome_dereplication"]["score"]["contamination"],
        #not in table
        N50_weight=config["genome_dereplication"]["score"]["N50"],
        size_weight=config["genome_dereplication"]["score"]["length"],
        opt_parameters=config["genome_dereplication"]["opt_parameters"],
        work_directory=lambda wc, output: os.path.dirname(output[0]),
    shell:
        " rm -rf {params.work_directory} ;"
        " dRep dereplicate "
        " --genomes {input[0]}/*.fasta "
        " --genomeInfo {input.quality} "
        " {params.filter} "
        " --length {params.filter_length} "
        " --completeness {params.filter_completeness} "
        " --contamination {params.filter_contamination} "
        " --S_ani {params.ANI} "
        " --cov_thresh {params.overlap} "
        " --completeness_weight {params.completeness_weight} "
        " --contamination_weight {params.contamination_weight} "
        " --N50_weight {params.N50_weight} "
        " --size_weight {params.size_weight} "
        " --processors {threads} "
        " {params.opt_parameters} "
        " {params.work_directory} "
        " &> {log} "


localrules:
    rename_genomes,


checkpoint rename_genomes:
    input:
        genomes="genomes/Dereplication/dereplicated_genomes",
    output:
        dir=directory("genomes/genomes"),
        mapfile_contigs="genomes/clustering/contig2genome.tsv",
        mapfile_genomes="genomes/clustering/old2newID.tsv",
        mapfile_bins="genomes/clustering/allbins2genome.tsv",
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


rule run_all_checkm_lineage_wf:
    input:
        touched_output="logs/checkm_init.txt",
        dir=genome_dir,
    output:
        "genomes/checkm/completeness.tsv",
        "genomes/checkm/storage/tree/concatenated.fasta",
    params:
        output_dir=lambda wc, output: os.path.dirname(output[0]),
        tmpdir=config["tmpdir"],
    conda:
        "%s/checkm.yaml" % CONDAENV
    threads: config.get("threads", 1)
    resources:
        time=config["runtime"]["long"],
        mem=config["large_mem"],
    log:
        "logs/genomes/checkm.log",
    benchmark:
        "logs/benchmarks/checkm_lineage_wf/all_genomes.tsv"
    shell:
        """
        rm -r {params.output_dir}
        checkm lineage_wf \
            --tmpdir {params.tmpdir} \
            --file {params.output_dir}/completeness.tsv \
            --tab_table \
            --quiet \
            --extension fasta \
            --threads {threads} \
            {input.dir} \
            {params.output_dir} &> {log}
        """


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
