import logging
import json
import pandas as pd
import os
import sys
from collections import defaultdict

from atlas import TAX_LEVELS
from atlas.utils import touch
from atlas.parsers import read_fasta


PROKKA_TSV_HEADER = ["contig_id", "gene_id", "ftype", "gene", "EC_number",
                     "product"]
REFSEQ_TSV_HEADER = ["contig", "orf", "taxonomy", "erfc", "orf_taxonomy",
                     "refseq_product", "refseq_evalue", "refseq_bitscore"]

EGGNOG_HEADER = ['query_name',
 'seed_eggNOG_ortholog',
 'seed_ortholog_evalue',
 'seed_ortholog_score',
 'predicted_gene_name',
 'GO_terms',
 'KEGG_KO',
 'BiGG_Reactions',
 'Annotation_tax_scope',
 'Matching_OGs',
 #'best_OG|evalue|score',
 'categories',
 #'eggNOG_HMM_model_annotation'
 ]


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
# #best_OG|evalue|score: Best matching Orthologous Groups (only in HMM mode)
# COG functional categories: COG functional category inferred from best matching OG
# eggNOG_HMM_model_annotation: eggNOG functional description inferred from best matching OG


# minus the sam/bam file name which is the last column
COUNTS_HEADER = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
MERGED_HEADER = ["contig_id", "gene_id", "ftype", "Length", "gene",
                 "EC_number", "product", "orf_taxonomy", "taxonomy", "erfc",
                 "count",
                 'seed_eggNOG_ortholog',
                 'seed_ortholog_score',
                 'predicted_gene_name',
                 'GO_terms',
                 'KEGG_KO',
                 'BiGG_Reactions',
                 'Matching_OGs',
                 'categories']
ANNOTATE_HEADER = ["contig_id", "gene_id", "ftype", "gene", "EC_number",
                   "product", "orf_taxonomy", "taxonomy", "erfc"]


BINNED_HEADER = MERGED_HEADER + ["bin_id", "checkm_bin_taxonomy_contained",
                 "checkm_bin_taxonomy_sister_lineage",
                 "checkm_bin_number_unique_markers", "checkm_bin_completeness",
                 "checkm_bin_contamination"]



def get_valid_dataframe(file_path, expected_cols, **kwargs):
    """
    Args:
        file_path (str): file path to data
        expected_cols (list): column headers necessary for validation

    Returns:
        pandas.core.frame.DataFrame

    Raises:
        ValueError: lists missing required columns
    """
    try:
        df = pd.read_csv(file_path, **kwargs)
    except UnicodeDecodeError:
        # UnicodeDecodeError: 'utf-8' codec can't decode byte 0x91 in position 35: invalid start byte
        kwargs["encoding"] = "iso-8859-1"
        df = pd.read_csv(file_path, **kwargs)
    missing_cols = []
    for c in expected_cols:
        if c not in df.columns:
            missing_cols.append(c)
    if len(missing_cols) > 0:
        raise ValueError("%s missing required columns: %s" % (file_path, ", ".join(missing_cols)))
    return df

## this is duplicated
# def do_merge(prokka_tsv, refseq_tsv, counts_tsv=None):
#     """Reads input files, creates temporary dataframes, and performs the merge."""
#     logging.info("Parsing Prokka TSV: %s" % prokka_tsv)
#     prokka_df = get_valid_dataframe(prokka_tsv, PROKKA_TSV_HEADER, sep="\t")
#
#     if counts_tsv:
#         logging.info("Parsing Counts TSV: %s" % counts_tsv)
#         counts_df = get_valid_dataframe(counts_tsv, COUNTS_HEADER, sep="\t",
#                         comment="#")
#         counts_df.rename(columns={counts_df.columns.tolist()[-1]:"count"},
#             inplace=True)
#
#         # merge prokka and counts
#         merged = pd.merge(left=counts_df, right=prokka_df, how="left",
#                      left_on=COUNTS_HEADER[0], right_on=PROKKA_TSV_HEADER[1])
#
#     logging.info("Parsing RefSeq file: %s" % refseq_tsv)
#     refseq_df = get_valid_dataframe(refseq_tsv, REFSEQ_TSV_HEADER, sep="\t")
#
#     if counts_tsv:
#         # merge in refseq data
#         merged = pd.merge(left=merged, right=refseq_df, how="left",
#                      left_on=COUNTS_HEADER[0], right_on=REFSEQ_TSV_HEADER[1])
#     else:
#         merged = pd.merge(left=prokka_df, right=refseq_df, how="left",
#                      left_on=PROKKA_TSV_HEADER[1],
#                      right_on=REFSEQ_TSV_HEADER[1])
#     return merged


def do_merge(prokka_tsv, refseq_tsv, counts_tsv=None,eggNOG=None):
    """Reads input files, creates temporary dataframes, and performs the merge."""
    logging.info("Parsing Prokka TSV: %s" % prokka_tsv)
    prokka_df = get_valid_dataframe(prokka_tsv, PROKKA_TSV_HEADER, sep="\t")

    if counts_tsv:
        logging.info("Parsing Counts TSV: %s" % counts_tsv)
        counts_df = get_valid_dataframe(counts_tsv, COUNTS_HEADER, sep="\t",
                        comment="#")
        counts_df.rename(columns={counts_df.columns.tolist()[-1]:"count"},
            inplace=True)

        # merge prokka and counts
        merged = pd.merge(left=counts_df, right=merged, how="left",
                     left_on=COUNTS_HEADER[0], right_on=PROKKA_TSV_HEADER[1])

    if eggNOG:
        logging.info("Parsing eggNOG TSV: %s" % counts_tsv)
        eggNOG_df = get_valid_dataframe(eggNOG, EGGNOG_HEADER, sep="\t")


        # merge prokka and counts
        merged = pd.merge(left=eggNOG_df, right=merged, how="left",
                     left_on='query_name', right_on=PROKKA_TSV_HEADER[1])


    logging.info("Parsing RefSeq file: %s" % refseq_tsv)
    refseq_df = get_valid_dataframe(refseq_tsv, REFSEQ_TSV_HEADER, sep="\t")

    if counts_tsv:
        # merge in refseq data
        merged = pd.merge(left=merged, right=refseq_df, how="left",
                     left_on=COUNTS_HEADER[0], right_on=REFSEQ_TSV_HEADER[1])
    else:
        merged = pd.merge(left=prokka_df, right=refseq_df, how="left",
                     left_on=PROKKA_TSV_HEADER[1],
                     right_on=REFSEQ_TSV_HEADER[1])
    return merged


def merge_bin_data(df, completeness_tsv, taxonomy_tsv, fasta_files):
    tax_df = pd.read_table(taxonomy_tsv)
    com_df = pd.read_table(completeness_tsv)
    m = pd.merge(left=tax_df, right=com_df, how="left", left_on="Bin Id",
                 right_on="Bin Id")

    contig_by_bin = defaultdict(list)
    for fasta in fasta_files:
        bin_id = os.path.splitext(os.path.basename(fasta))[0]
        with open(fasta) as fh:
            for name, seq in read_fasta(fh):
                contig_by_bin[bin_id].append(name)

    contig_dict = {}
    for bin_id, contig_ids in contig_by_bin.items():
        for contig_id in contig_ids:
            contig_dict[contig_id] = {"bin_id":bin_id,
                "checkm_bin_taxonomy_contained":m.loc[m["Bin Id"] == bin_id]["Taxonomy (contained)"].values[0],
                "checkm_bin_taxonomy_sister_lineage":m.loc[m["Bin Id"] == bin_id]["Taxonomy (sister lineage)"].values[0],
                "checkm_bin_number_unique_markers":m.loc[m["Bin Id"] == bin_id]["# unique markers (of 43)"].values[0],
                "checkm_bin_completeness":m.loc[m["Bin Id"] == bin_id]["Completeness"].values[0],
                "checkm_bin_contamination":m.loc[m["Bin Id"] == bin_id]["Contamination"].values[0]
            }

    contig_df = pd.DataFrame.from_dict(contig_dict, orient="index")
    contig_df = contig_df.reset_index().rename(columns={"index":"binned_contig_id"})

    merged = pd.merge(left=df, right=contig_df, how="left", left_on="Chr",
                      right_on="binned_contig_id")
    return merged


def merge_tables(prokka_tsv, refseq_tsv, output,counts_tsv=None, eggNOG= None,
    completeness=None, taxonomy=None, fastas=None):
    """

    Count data is a TSV formatted with a header:

        \b
        Geneid Chr Start  End Strand Length /path/example.bam
        orf1     1     1  500      +    500                50
        orf2     1   601  900      +    300               300
        orf3     1  1201 1500      +    300               200
    """
    if counts_tsv and eggNOG and  completeness and taxonomy and fastas:
        df = do_merge(prokka_tsv, refseq_tsv, counts_tsv,eggNOG)
        df = merge_bin_data(df, completeness, taxonomy, fastas)
        df.to_csv(output, sep="\t", columns=BINNED_HEADER, index=False)
    else:
        df = do_merge(prokka_tsv, refseq_tsv, counts_tsv, eggNOG)
        if counts_tsv:
            df.to_csv(output, sep="\t", columns=MERGED_HEADER, index=False)
        else:
            df.to_csv(output, sep="\t", columns=ANNOTATE_HEADER, index=False)
