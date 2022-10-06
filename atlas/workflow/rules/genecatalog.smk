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
        "Genecatalog/counts/sequence_infos.tsv",
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
        "Genecatalog/gene_catalog.fna",
    output:
        "ref/Genecatalog.mmi",
    log:
        "logs/Genecatalog/alignment/index.log",
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


rule align_reads_to_Genecatalog:
    input:
        mmi=rules.index_genecatalog.output,
        reads=lambda wc: get_quality_controlled_reads(wc, include_se=True),
    output:
        temp("Genecatalog/alignments/{sample}.bam"),
    log:
        "logs/Genecatalog/alignment/{sample}_map.log",
    threads: config["threads"]
    resources:
        mem=config["mem"],
        mem_mb=config["mem"] * 1000,
    conda:
        "../envs/minimap.yaml"
    shell:
        " (cat {input.reads} | "
        " minimap2 "
        " -t {threads} "
        "-a "
        "-x sr "
        " {input.mmi} "
        " - "
        " | samtools view -b "
        " > {output} ) 2>{log}"


rule pileup_Genecatalog:
    input:
        bam=rules.align_reads_to_Genecatalog.output,
    output:
        covstats=temp("Genecatalog/alignments/{sample}_coverage.tsv"),
        rpkm=temp("Genecatalog/alignments/{sample}_rpkm.tsv"),
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
        " -Xmx{resources.java_mem}G "
        " 2> {log} "


rule combine_gene_coverages:
    input:
        covstats=expand("Genecatalog/alignments/{sample}_coverage.tsv", sample=SAMPLES),
    output:
        "Genecatalog/counts/median_coverage.parquet",
        "Genecatalog/counts/Nmapped_reads.parquet",
    log:
        "logs/Genecatalog/counts/combine_gene_coverages.log",
    params:
        samples=SAMPLES,
    threads: 1
    resources:
        mem=config["large_mem"],
    script:
        "../scripts/combine_gene_coverages.py"


# TODO: combine RPKM

# TODO: caluclate mapping rate from pileup mapping files
# logs/Genecatalog/alignment/sample2_pileup.log
# Reads:                                  1207217
# Mapped reads:                           1071071


###########
## EGG NOG
##########

# # this rule specifies the more general eggNOG rules

# output with wildcards "{folder}/{prefix}.emapper.tsv"


rule eggNOG_homology_search:
    input:
        eggnog_db_files=get_eggnog_db_file(),
        faa="{folder}/{prefix}.faa",
    output:
        temp("{folder}/{prefix}.emapper.seed_orthologs"),
    params:
        data_dir=EGGNOG_DIR,
        prefix="{folder}/{prefix}",
    resources:
        mem=config["mem"],
    threads: config["threads"]
    shadow:
        "minimal"
    conda:
        "%s/eggNOG.yaml" % CONDAENV
    log:
        "{folder}/logs/{prefix}/eggNOG_homology_search_diamond.log",
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
        temp("{folder}/{prefix}.emapper.annotations"),
    params:
        data_dir=(
            config["virtual_disk"] if config["eggNOG_use_virtual_disk"] else EGGNOG_DIR
        ),
        prefix="{folder}/{prefix}",
        copyto_shm="t" if config["eggNOG_use_virtual_disk"] else "f",
    threads: config.get("threads", 1)
    resources:
        mem=calculate_mem_eggnog(),
    shadow:
        "minimal"
    conda:
        "%s/eggNOG.yaml" % CONDAENV
    log:
        "{folder}/logs/{prefix}/eggNOG_annotate_hits_table.log",
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


#
# localrules: get_Genecatalog_annotations
# rule get_Genecatalog_annotations:
#     input:
#         Genecatalog= 'Genecatalog/gene_catalog.fna".fna',
#         eggNOG= expand('{sample}/annotation/eggNOG.tsv',sample=SAMPLES),
#         refseq= expand('{sample}/annotation/refseq/{sample}_tax_assignments.tsv',sample=SAMPLES),
#         scg= expand("Genecatalog/annotation/single_copy_genes_{domain}.tsv",domain=['bacteria','archaea'])
#     output:
#         annotations= "Genecatalog/annotations.tsv",
#     run:
#         import pandas as pd
#
#         gene_ids=[]
#         with open(input.Genecatalog) as fasta_file:
#             for line in fasta_file:
#                 if line[0]=='>':
#                     gene_ids.append(line[1:].strip().split()[0])
#
#         eggNOG=pd.DataFrame()
#         for annotation_file in input.eggNOG:
#             eggNOG=eggNOG.append(pd.read_csv(annotation_file, index_col=0,sep='\t'))
#
#         refseq=pd.DataFrame()
#         for annotation_file in input.refseq:
#             refseq=refseq.append(pd.read_csv(annotation_file, index_col=1,sep='\t'))
#
#         scg=pd.DataFrame()
#         for annotation_file in input.scg:
#             d= pd.read_csv(annotation_file, index_col=0,header=None,sep='\t')
#             d.columns = 'scg_'+ os.path.splitext(annotation_file)[0].split('_')[-1] # bacteria or archaea
#             scg=scg.append(d)
#
#
#         annotations= refseq.join(eggNOG).join(scg).loc[gene_ids]
#         annotations.to_csv(output.annotations,sep='\t')


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


localrules:
    gene_subsets,
    combine_egg_nogg_annotations,


checkpoint gene_subsets:
    input:
        "Genecatalog/gene_catalog.faa",
    output:
        directory("Genecatalog/subsets/genes"),
    params:
        subset_size=config["genecatalog"]["SubsetSize"],
    conda:
        "../envs/sequence_utils.yaml"
    log:
        "logs/Genecatalog/clustering/split_genecatalog.log",
    script:
        "../scripts/split_genecatalog.py"


def combine_genecatalog_annotations_input(wildcards):
    dir_for_subsets = checkpoints.gene_subsets.get(**wildcards).output[0]
    (Subset_names,) = glob_wildcards(os.path.join(dir_for_subsets, "{subset}.faa"))
    return expand(
        "Genecatalog/subsets/genes/{subset}.emapper.annotations", subset=Subset_names
    )


rule combine_egg_nogg_annotations:
    input:
        combine_genecatalog_annotations_input,
    output:
        "Genecatalog/annotations/eggNog.tsv.gz",
    run:
        import pandas as pd

        # read input files one after the other
        for i, annotation_table in enumerate(input):
            D = pd.read_csv(annotation_table, header=None, sep="\t")
            # Add headers, to verify size
            D.columns = EGGNOG_HEADER
            # appedn to output file, header only the first time
            D.to_csv(
                output[0],
                sep="\t",
                index=False,
                header=(i == 0),
                compression="gzip",
                mode="a",
            )


rule gene2genome:
    input:
        contigs2bins="genomes/clustering/all_contigs2bins.tsv.gz",
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
