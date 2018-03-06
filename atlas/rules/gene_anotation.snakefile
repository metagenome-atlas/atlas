

faa = "{sample}/annotation/{sample}_genes.faa"

CONDAENV = "../envs"
#eggNOG_DATABASES =  ["none"] #bact arch viruses none for only diamond " ".join(config.get("eggNOG_databases", eggNOG_DATABASES))

EGGNOG_DIR  = os.path.join(os.path.dirname(os.path.realpath(config.get("diamond_db", "."))), "eggNOG")


rule all:
    input:
        "combined/annotation/eggNOG/combined.emapper.annotations"


# rule gene_calling:
#     input:
#         "{folder}/{sample}_contigs.fasta"
#     output:
#         fna="{folder}/annotation/{sample}_genes.fna",
#         faa="{folder}/annotation/{sample}_genes.faa",
#         gff="{folder}/annotation/{sample}_genes.gff"
#
#     conda:
#         "%s/gene_catalog.yaml" % CONDAENV
#     log:
#         "{folder}/logs/predict_genes_{sample}.log"
#     threads:
#         1
#     shell:
#         """
#             prodigal -i {input} -o {output.gff} -d {output.fna} -a {output.faa} -p meta -f gff 2> >(tee {log})
#         """





#eggNOG_DATABASES = expand("{data_dir}/{database_files}",data_dir= EGGNOG_DIR,database_files= ["eggnog_proteins.dmnd","eggnog.db","og2level.tsv","OG.fata"
localrules: download_eggnog_data

rule download_eggnog_data:
    output:
        touch("%s/download_eggnog_data.sucess" % EGGNOG_DIR)
    params:
        dbs = "none",
        data_dir = EGGNOG_DIR
    threads:
        2
    conda:
        "%s/eggNOG.yaml" % CONDAENV
    log:
        "logs/download/download_eggnog_data.log"
    shell:
        """
            download_eggnog_data.py -y -f -q --data_dir {params.data_dir} {params.dbs}
        """

#HIGHT troughput : split faa in 1Mio faa chunks for next step
localrules: eggNOG_homology_search
rule eggNOG_homology_search:
    input:
        "%s/download_eggnog_data.sucess" % EGGNOG_DIR ,
        faa = faa,
    output:
        "{sample}/annotation/eggNOG/{sample}.emapper.seed_orthologs",
    params:
        data_dir = EGGNOG_DIR,
        prefix = lambda wc,output: output[0].replace(".emapper.seed_orthologs","")
    threads:
        config.get("threads", 1)
    conda:
        "%s/eggNOG.yaml" % CONDAENV
    log:
        "{sample}/logs/eggNOG/homology_search.log"
    shell:
        """
            emapper.py -m diamond --no_annot --no_file_comments --data_dir {params.data_dir} \
            --cpu {threads} -i {input.faa} -o {params.prefix} --override 2> >(tee {log})
        """


#HIGHT troughput : concat emapper.seed_orthologs chunks
# run on single machine
rule eggNOG_annotation:
    input:
        "%s/download_eggnog_data.sucess" % EGGNOG_DIR ,
        seed = rules.eggNOG_homology_search.output
    output:
        "{sample}/annotation/eggNOG/{sample}.emapper.annotations"
    params:
        data_dir = EGGNOG_DIR,
        prefix = lambda wc,output: output[0].replace(".emapper.annotations","")
    threads:
        config.get("threads", 1)
    conda:
        "%s/eggNOG.yaml" % CONDAENV
    log:
        "{sample}/logs/eggNOG/annotation.log"
    shell:
        """
            emapper.py --annotate_hits_table {input.seed} --no_file_comments \
            --override -o {params.prefix} --cpu {threads} --data_dir {params.data_dir} 2> >(tee {log})
        """


# M.columns ['query_name',
#  'seed_eggNOG_ortholog',
#  'seed_ortholog_evalue',
#  'seed_ortholog_score',
#  'predicted_gene_name',
#  'GO_terms',
#  'KEGG_KO',
#  'BiGG_Reactions',
#  'Annotation_tax_scope',
#  'Matching_OGs',
#  'best_OG|evalue|score',
#  'categories',
#  'eggNOG_HMM_model_annotation']

# This file provides final annotations of each query. Tab-delimited columns in the file are:

# query_name: query sequence name
# seed_eggNOG_ortholog: best protein match in eggNOG
# seed_ortholog_evalue: best protein match (e-value)
# seed_ortholog_score: best protein match (bit-score)
# predicted_gene_name: Predicted gene name for query sequences
# GO_terms: Comma delimited list of predicted Gene Ontology terms
# KEGG_KO: Comma delimited list of predicted KEGG KOs
# BiGG_Reactions: Comma delimited list of predicted BiGG metabolic reactions
# Annotation_tax_scope: The taxonomic scope used to annotate this query sequence
# Matching_OGs: Comma delimited list of matching eggNOG Orthologous Groups
# best_OG|evalue|score: Best matching Orthologous Groups (only in HMM mode)
# COG functional categories: COG functional category inferred from best matching OG
# eggNOG_HMM_model_annotation: eggNOG functional description inferred from best matching OG
