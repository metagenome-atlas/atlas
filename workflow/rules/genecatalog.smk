import os


rule filter_genes:
    input:
        fna="{gene_file}.fna",
        faa="{gene_file}.faa",
    output:
        fna=temp("{gene_file}.filtered.fna"),
        faa=temp("{gene_file}.filtered.faa"),
        short=temp("{gene_file}.short.faa"),
    conda:
        "../envs/fasta.yaml"
    threads: 1
    log:
        "logs/Genecatalog/filter_genes/{gene_file}.log",
    params:
        minlength_nt=config["genecatalog"]["minlength_nt"],
    script:
        "../scripts/filter_genes.py"


if config["genecatalog"]["source"] == "contigs":

    localrules:
        concat_genes,

    rule concat_genes:
        input:
            faa=expand(
                "{sample}/annotation/predicted_genes/{sample}.filtered.faa",
                sample=SAMPLES,
            ),
            fna=expand(
                "{sample}/annotation/predicted_genes/{sample}.filtered.fna",
                sample=SAMPLES,
            ),
            short=expand(
                "{sample}/annotation/predicted_genes/{sample}.short.faa",
                sample=SAMPLES,
            ),
        output:
            faa=temp("Genecatalog/all_genes/predicted_genes.faa"),
            fna=temp("Genecatalog/all_genes/predicted_genes.fna"),
            short=temp("Genecatalog/all_genes/short_genes.faa"),
        run:
            from utils.io import cat_files

            cat_files(input.faa, output.faa)
            cat_files(input.fna, output.fna)
            cat_files(input.short, output.short)


else:

    localrules:
        concat_genes,

    rule concat_genes:
        input:
            "genomes/annotations/orf2genome.tsv",
            faa=lambda wc: get_all_genes(wc, ".filtered.faa"),
            fna=lambda wc: get_all_genes(wc, ".filtered.fna"),
            short=lambda wc: get_all_genes(wc, ".short.faa"),
        output:
            faa=temp("Genecatalog/all_genes/predicted_genes.faa"),
            fna=temp("Genecatalog/all_genes/predicted_genes.fna"),
            short=temp("Genecatalog/all_genes/short_genes.faa"),
        run:
            from utils.io import cat_files

            cat_files(input.faa, output.faa)
            cat_files(input.fna, output.fna)
            cat_files(input.short, output.short)


if (config["genecatalog"]["clustermethod"] == "linclust") or (
    config["genecatalog"]["clustermethod"] == "mmseqs"
):

    rule cluster_genes:
        input:
            faa="Genecatalog/all_genes/predicted_genes.faa",
        output:
            db=temp(directory("Genecatalog/all_genes/predicted_genes")),
            clusterdb=temp(directory("Genecatalog/clustering/mmseqs")),
        conda:
            "%s/mmseqs.yaml" % CONDAENV
        log:
            "logs/Genecatalog/clustering/cluster_proteins.log",
        threads: config.get("threads", 1)
        params:
            tmpdir=os.path.join(config["tmpdir"], "mmseqs"),
            clustermethod=(
                "linclust"
                if config["genecatalog"]["clustermethod"] == "linclust"
                else "cluster"
            ),
            coverage=config["genecatalog"]["coverage"],  #0.8,
            minid=config["genecatalog"]["minid"],  # 0.00
            extra=config["genecatalog"]["extra"],
            clusterdb=lambda wc, output: os.path.join(output.clusterdb, "clusterdb"),
            db=lambda wc, output: os.path.join(output.db, "inputdb"),
        shell:
            """
            mkdir -p {params.tmpdir} {output} 2>> {log}
            mmseqs createdb {input.faa} {params.db} &> {log}

            mmseqs {params.clustermethod} -c {params.coverage} \
            --min-seq-id {params.minid} {params.extra} \
            --threads {threads} {params.db} {params.clusterdb} {params.tmpdir}  &>>  {log}

            rm -fr  {params.tmpdir} 2>> {log}
            """

    rule get_rep_proteins:
        input:
            db=rules.cluster_genes.output.db,
            clusterdb=rules.cluster_genes.output.clusterdb,
        output:
            cluster_attribution=temp("Genecatalog/orf2gene_oldnames.tsv"),
            rep_seqs_db=temp(directory("Genecatalog/protein_catalog")),
            rep_seqs=temp("Genecatalog/representatives_of_clusters.faa"),
        conda:
            "%s/mmseqs.yaml" % CONDAENV
        log:
            "logs/Genecatalog/clustering/get_rep_proteins.log",
        threads: config.get("threads", 1)
        params:
            clusterdb=lambda wc, input: os.path.join(input.clusterdb, "clusterdb"),
            db=lambda wc, input: os.path.join(input.db, "inputdb"),
        shell:
            """
            mmseqs createtsv {params.db} {params.db} {params.clusterdb} {output.cluster_attribution}  &> {log}

            mkdir {output.rep_seqs_db} 2>> {log}

            mmseqs result2repseq {params.db} {params.clusterdb} {output.rep_seqs_db}/db  &>> {log}

            mmseqs result2flat {params.db} {params.db} {output.rep_seqs_db}/db {output.rep_seqs}  &>> {log}

            """

    rule get_cds_of_proteins:
        input:
            all="Genecatalog/all_genes/predicted_genes.fna",
            names="Genecatalog/representatives_of_clusters.faa",
        output:
            temp("Genecatalog/representatives_of_clusters.fna"),
        conda:
            "../envs/required_packages.yaml"
        threads: 1
        resources:
            mem=config["mem"],
            java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
        log:
            "logs/Genecatalog/clustering/get_cds_of_proteins.log",
        shell:
            " filterbyname.sh "
            " in={input.all}"
            " names={input.names}"
            " include=t"
            " out={output} "
            " -Xmx{resources.java_mem}G "
            " 2> {log}"

    rule generate_orf_info:
        input:
            cluster_attribution="Genecatalog/orf2gene_oldnames.tsv",
        output:
            cluster_attribution="Genecatalog/clustering/orf_info.parquet",
            rep2genenr="Genecatalog/clustering/representative2genenr.tsv",
        threads: 1
        log:
            "logs/Genecatalog/clustering/generate_orf_info.log",
        script:
            "../scripts/generate_orf_info.py"
# cluster genes with cd-hit-est



elif config["genecatalog"]["clustermethod"] == "cd-hit-est":

    include: "cdhit.smk"


else:
    raise Exception(
        "Didn't understood the genecatalog clustermethod: {}".format(
            config["genecatalog"]["clustermethod"]
        )
    )


localrules:
    rename_gene_catalog,


rule rename_gene_catalog:
    input:
        fasta="Genecatalog/representatives_of_clusters.{ext}",
        rep2genenr="Genecatalog/clustering/representative2genenr.tsv",
    output:
        "Genecatalog/gene_catalog.{ext}",
    log:
        "logs/Genecatalog/clustering/rename_gene_catalog_{ext}.log",
    script:
        "../scripts/rename_genecatalog.py"


rule get_genecatalog_seq_info:
    input:
        "Genecatalog/gene_catalog.fna",
    output:
        temp("Genecatalog/counts/sequence_infos.tsv"),
    log:
        "logs/Genecatalog/get_seq_info.log",
    conda:
        "../envs/required_packages.yaml"
    threads: 1
    resources:
        mem=config["simplejob_mem"],
        java_mem=int(config["simplejob_mem"] * JAVA_MEM_FRACTION),
    shell:
        "stats.sh gcformat=4 gc={output} in={input} &> {log}"


rule index_genecatalog:
    input:
        target="Genecatalog/gene_catalog.fna",
    output:
        temp("ref/Genecatalog.mmi"),
    log:
        "logs/Genecatalog/alignment/index.log",
    params:
        index_size="12G",
    wrapper:
        "v1.19.0/bio/minimap2/index"


rule concat_all_reads:
    input:
        lambda wc: get_quality_controlled_reads(wc, include_se=True),
    output:
        temp("Intermediate/genecatalog/alignments/{sample}.fastq.gz")
    log:
        "logs/Genecatalog/alignment/concat_reads/{sample}.log"
    threads:
        1
    resources:
        mem_mb=300
    shell:
        "cat {input} > {output} 2> {log}"
    

rule align_reads_to_Genecatalog:
    input:
        target=rules.index_genecatalog.output,
        query= rules.concat_all_reads.output[0]
    output:
        temp("Genecatalog/alignments/{sample}.bam"),
    log:
        "logs/Genecatalog/alignment/{sample}_map.log",
    threads: config["threads"]
    resources:
        mem_mb=config["mem"] * 1000,
    params:
        extra="-x sr --split-prefix {sample}_split_ ",
        sort="coordinate",
    wrapper:
        "v1.19.0/bio/minimap2/aligner"



rule pileup_Genecatalog:
    input:
        bam=rules.align_reads_to_Genecatalog.output,
    output:
        covstats=temp("Genecatalog/alignments/{sample}_coverage.tsv"),
        rpkm=temp("Genecatalog/alignments/{sample}_rpkm.tsv"),
    params:
        minmapq=config["minimum_map_quality"]
    log:
        "logs/Genecatalog/alignment/{sample}_pileup.log",
    conda:
        "../envs/required_packages.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
    shell:
        " pileup.sh "
        " in={input.bam}"
        " covstats={output.covstats} "
        " rpkm={output.rpkm} "
        " secondary=t "
        " minmapq={params.minmapq} "
        " -Xmx{resources.java_mem}G "
        " 2> {log} "


rule gene_pileup_as_parquet:
    input:
        cov="Genecatalog/alignments/{sample}_coverage.tsv",
        #rpkm = "Genecatalog/alignments/{sample}_rpkm.tsv"
    output:
        "Genecatalog/alignments/{sample}_coverage.parquet",
    threads: 1
    resources:
        mem=config["simplejob_mem"],
        time_min=config["runtime"]["simplejob"]*60,
    log:
        "logs/Genecatalog/counts/parse_gene_coverages/{sample}.log",
    run:
        try:
            import pandas as pd
            from utils.parsers_bbmap import read_pileup_coverage

            data = read_pileup_coverage(
                input[0],
                coverage_measure="Median_fold",
                other_columns=["Avg_fold", "Covered_percent", "Read_GC", "Std_Dev"],
            )
            data.index.name = "GeneName"
            data.sort_index(inplace=True)

            # rpkm = pd.read_csv(input[1],sep='\t',skiprows=4,usecols=["#Name","RPKM"],index_col=0).sort_index()

            data.reset_index().to_parquet(output[0])

        except Exception as e:
            import traceback

            with open(log[0], "w") as logfile:
                traceback.print_exc(file=logfile)

            raise e


def get_combine_cov_time():
    estimated_time = (0.5 * len(SAMPLES) + 20) / 60

    config_time = config["runtime"]["long"]

    if config_time < estimated_time:
        logger.error(
            f"For rule combine_gene_coverages we estimate 0.5 min/ per sample = {estimated_time:.1f} h. "
            "You provided to little. \n Increase time in config file: \nruntime:\n  long\n."
        )
        raise Exception("Not long enough runtime provided. ")

    return config_time * 60


rule combine_gene_coverages:
    input:
        covstats=expand(
            "Genecatalog/alignments/{sample}_coverage.parquet", sample=SAMPLES
        ),
        info="Genecatalog/counts/sequence_infos.tsv",
    output:
        cov="Genecatalog/counts/median_coverage.h5",
        counts="Genecatalog/counts/Nmapped_reads.h5",
        sample_info="Genecatalog/counts/sample_coverage_stats.tsv",
        gene_info="Genecatalog/counts/gene_coverage_stats.parquet",
    log:
        "logs/Genecatalog/counts/combine_gene_coverages.log",
    params:
        samples=SAMPLES,
    conda:
        "../envs/hdf.yaml"
    threads: 1
    resources:
        mem=config["simplejob_mem"],
        time_min=get_combine_cov_time(),
    script:
        "../scripts/combine_gene_coverages.py"


# TODO: combine RPKM

# TODO: caluclate mapping rate from pileup mapping files
# logs/Genecatalog/alignment/sample2_pileup.log
# Reads:                                  1207217
# Mapped reads:                           1071071


# Temporary
# move old subset folder to new
old_subset_folder = Path("Genecatalog/subsets/genes")
new_subset_folder = "Intermediate/genecatalog/subsets"
if old_subset_folder.exists():

    logger.info(f"I move {old_subset_folder} to {new_subset_folder}")

    import shutil

    shutil.move(old_subset_folder, new_subset_folder)


localrules:
    gene_subsets,


checkpoint gene_subsets:
    input:
        "Genecatalog/gene_catalog.faa",
    output:
        directory("Intermediate/genecatalog/subsets"),
    params:
        subset_size=config["genecatalog"]["SubsetSize"],
    conda:
        "../envs/sequence_utils.yaml"
    log:
        "logs/Genecatalog/clustering/split_genecatalog.log",
    script:
        "../scripts/split_genecatalog.py"


def get_subset_names(wildcards):

    dir_for_subsets = Path(checkpoints.gene_subsets.get(**wildcards).output[0])
    subset_names = glob_wildcards(str(dir_for_subsets / "{subset}.faa")).subset

    return subset_names


###########
## EGG NOG
##########

# # this rule specifies the more general eggNOG rules

# output with wildcards "{folder}/{prefix}.emapper.tsv"


rule eggNOG_homology_search:
    input:
        eggnog_db_files=get_eggnog_db_file(),
        faa="Intermediate/genecatalog/subsets/{subset}.faa",
    output:
        temp(
            "Intermediate/genecatalog/annotation/eggNOG/{subset}.emapper.seed_orthologs"
        ),
    params:
        data_dir=EGGNOG_DIR,
        prefix=lambda wc, output: output[0].replace(".emapper.seed_orthologs", ""),
    resources:
        mem=config["mem"],
    threads: config["threads"]
    shadow:
        "minimal"
    conda:
        "../envs/eggNOG.yaml"
    log:
        "logs/genecatalog/annotation/eggnog/{subset}_homology_search_diamond.log",
    shell:
        """
        emapper.py -m diamond --no_annot --no_file_comments \
            --data_dir {params.data_dir} --cpu {threads} -i {input.faa} \
            -o {params.prefix} --override 2> {log}
        """


def calculate_mem_eggnog():
    return 2 * config["simplejob_mem"] + (
        37 if config["eggNOG_use_virtual_disk"] else 0
    )


rule eggNOG_annotation:
    input:
        eggnog_db_files=get_eggnog_db_file(),
        seed=rules.eggNOG_homology_search.output,
    output:
        temp("Intermediate/genecatalog/annotation/eggNOG/{subset}.emapper.annotations"),
    params:
        data_dir=(
            config["virtual_disk"] if config["eggNOG_use_virtual_disk"] else EGGNOG_DIR
        ),
        prefix=lambda wc, output: output[0].replace(".emapper.annotations", ""),
        copyto_shm="t" if config["eggNOG_use_virtual_disk"] else "f",
    threads: config.get("threads", 1)
    resources:
        mem=calculate_mem_eggnog(),
    shadow:
        "minimal"
    conda:
        "../envs/eggNOG.yaml"
    log:
        "logs/genecatalog/annotation/eggnog/{subset}_annotate_hits_table.log",
    shell:
        """

        if [ {params.copyto_shm} == "t" ] ;
        then
            cp {EGGNOG_DIR}/eggnog.db {params.data_dir}/eggnog.db 2> {log}
            cp {EGGNOG_DIR}/eggnog_proteins.dmnd {params.data_dir}/eggnog_proteins.dmnd 2>> {log}
        fi

        emapper.py --annotate_hits_table {input.seed} --no_file_comments \
          --override -o {params.prefix} --cpu {threads} --data_dir {params.data_dir} 2>> {log}
        """


def combine_eggnog_annotations(wildcards):
    all_subsets = get_subset_names(wildcards)

    return expand(rules.eggNOG_annotation.output[0], subset=all_subsets)


rule combine_egg_nogg_annotations:
    input:
        combine_eggnog_annotations,
    output:
        "Genecatalog/annotations/eggNOG.parquet",
    log:
        "logs/genecatalog/annotation/eggNOG/combine.log",
    resources:
        time=config["runtime"]["default"],
    run:
        try:

            import pandas as pd

            Tables = [
                pd.read_csv(file, index_col=None, header=None, sep="\t")
                for file in input
            ]

            combined = pd.concat(Tables, axis=0)

            del Tables

            combined.columns = EGGNOG_HEADER

            #           combined.sort_values("Gene",inplace=True)

            combined.to_parquet(output[0], index=False)
        except Exception as e:

            import traceback

            with open(log[0], "w") as logfile:
                traceback.print_exc(file=logfile)

            raise e


## Temporary
rule convert_eggNOG_tsv2parquet:
    input:
        "Genecatalog/annotations/eggNog.tsv.gz",
    output:
        "Genecatalog/annotations/eggNOG.parquet",
    resources:
        time=config["runtime"]["default"],
    log:
        "logs/genecatalog/annotation/eggNOG/tsv2parquet.log",
    run:
        try:
            import pandas as pd

            df = pd.read_table(input[0], index_col=None)

            df.to_parquet(output[0], index=False)

        except Exception as e:

            import traceback

            with open(log[0], "w") as logfile:
                traceback.print_exc(file=logfile)

            raise e


ruleorder: convert_eggNOG_tsv2parquet > combine_egg_nogg_annotations


##############################################
### DRAM
###############################################


rule DRAM_annotate_genecatalog:
    input:
        faa="Intermediate/genecatalog/subsets/{subset}.faa",
        config=get_dram_config,
    output:
        annotations=temp(
            "Intermediate/genecatalog/annotation/dram/{subset}/annotations.tsv"
        ),
        genes=temp("Intermediate/genecatalog/annotation/dram/{subset}/genes.faa"),
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["long"],
    conda:
        "../envs/dram.yaml"
    params:
        extra=config.get("dram_extra", ""),
        outdir=lambda wc, output: Path(output[0]).parent,
    log:
        "logs/Genecatalog/annotation/dram/{subset}.log",
        "logs/Genecatalog/annotation/dram/{subset}.logfile",
    shell:
        " rm -rf {params.outdir} &> {log[0]};"
        "\n"
        " DRAM.py annotate_genes "
        " --input_faa {input.faa}"
        " --config_loc {input.config} "
        " --output_dir {params.outdir} "
        " --threads {threads} "
        " {params.extra} "
        " --log_file_path {log[1]} "
        " --verbose &>> {log[0]}"


def combine_genecatalog_dram_input(wildcards):

    all_subsets = get_subset_names(wildcards)

    return expand(
        rules.DRAM_annotate_genecatalog.output.annotations, subset=all_subsets
    )


rule combine_dram_genecatalog_annotations:
    input:
        combine_genecatalog_dram_input,
    output:
        directory("Genecatalog/annotations/dram"),
    resources:
        time=config["runtime"]["default"],
    log:
        "logs/genecatalog/annotation/dram/combine.log",
    script:
        "../scripts/combine_dram_gene_annotations.py"


rule gene2genome:
    input:
        contigs2bins="Binning/{binner}/contigs2bins.tsv.gz".format(
            binner=config["final_binner"]
        ),
        contigs2mags="genomes/clustering/contig2genome.tsv",
        old2newID="genomes/clustering/old2newID.tsv",
        orf_info="Genecatalog/clustering/orf_info.parquet",
    params:
        renamed_contigs=config["rename_mags_contigs"]
        & (config["genecatalog"]["source"] == "contigs"),
    output:
        "genomes/annotations/gene2genome.parquet",
    log:
        "logs/genomes/annotations/gene2genome.log",
    script:
        "../scripts/gene2genome.py"


# after combination need to add eggNOG headerself.
# "{folder}/{prefix}_eggNOG.tsv"
#
# ############## Canopy clustering
#
# rule reformat_for_canopy:
#         input:
#             "mapresults/Genecatalog_CE/combined_Nmaped_reads.tsv"
#         output:
#             "mapresults/Genecatalog_CE/nseq.tsv"
#         run:
#             import pandas as pd
#
#             D= pd.read_csv(input[0], index_col=0,sep='\t')
#             D.index= D.index.map(lambda s: s.split()[0])
#             D=D.astype(int)
#             D.to_csv(output[0],sep='\t',header=False)
#
#
# rule canopy_clustering:
#     input:
#         rules.reformat_for_canopy.output
#     output:
#         cluster="mapresults/Genecatalog_CE/canopy_cluster.tsv",
#         profile="mapresults/Genecatalog_CE/cluster_profiles.tsv",
#     params:
#         canopy_params=config.get("canopy_params","")
#     log:
#         "mapresults/Genecatalog_CE/canopy.log"
#     benchmark:
#         "logs/benchmarks/canopy_clustering.txt"
#     conda:
#         "%s/canopy.yaml" % CONDAENV
#     threads:
#         12
#     resources:
#         mem= 220
#     shell:
#         """
#         canopy -i {input} -o {output.cluster} -c {output.profile} -n {threads} --canopy_size_stats_file {log} {params.canopy_params} 2> {log}
#
#         """
rule predict_single_copy_genes:
    input:
        "Genecatalog/gene_catalog.faa",
    output:
        "Genecatalog/annotation/single_copy_genes_{domain}.tsv",
    params:
        script_dir=os.path.dirname(os.path.abspath(workflow.snakefile)),
        key=lambda wc: wc.domain[:3],  #bac for bacteria, #arc for archaea
    conda:
        "%s/DASTool.yaml" % CONDAENV  # needs pearl
    log:
        "logs/Genecatalog/annotation/predict_single_copy_genes_{domain}.log",
    shadow:
        "shallow"
    threads: config["threads"]
    shell:
        " DIR=$(dirname $(readlink -f $(which DAS_Tool))) "
        ";"
        " ruby {params.script_dir}/rules/scg_blank_diamond.rb diamond"
        " {input} "
        " $DIR\/db/{params.key}.all.faa "
        " $DIR\/db/{params.key}.scg.faa "
        " $DIR\/db/{params.key}.scg.lookup "
        " {threads} "
        " 2> {log} "
        " ; "
        " mv {input[0]}.scg {output}"
