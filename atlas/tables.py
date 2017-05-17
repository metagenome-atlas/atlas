import logging
import json
import pandas as pd
import sys
from atlas import TAX_LEVELS
from atlas.utils import touch


PROKKA_TSV_HEADER = ["contig_id", "locus_tag", "ftype", "gene", "EC_number", "product"]
REFSEQ_TSV_HEADER = ["contig", "orf", "taxonomy", "erfc", "orf_taxonomy", "refseq_product",
                     "refseq_evalue", "refseq_bitscore"]
# minus the sam/bam file name which is the last column
COUNTS_HEADER = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
MERGED_HEADER = ["contig_id", "locus_tag", "ftype", "Length", "gene", "EC_number", "product",
                 "orf_taxonomy", "taxonomy", "erfc", "count"]
ANNOTATE_HEADER = ["contig_id", "locus_tag", "ftype", "gene", "EC_number", "product",
                   "orf_taxonomy", "taxonomy", "erfc"]


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


def do_merge(prokka_tsv, refseq_tsv, counts_tsv=None):
    """Reads input files, creates temporary dataframes, and performs the merge."""
    logging.info("Parsing Prokka TSV: %s" % prokka_tsv)
    prokka_df = get_valid_dataframe(prokka_tsv, PROKKA_TSV_HEADER, sep="\t")

    if counts_tsv:
        logging.info("Parsing Counts TSV: %s" % counts_tsv)
        counts_df = get_valid_dataframe(counts_tsv, COUNTS_HEADER, sep="\t", comment="#")
        counts_df.rename(columns={counts_df.columns.tolist()[-1]:"count"}, inplace=True)

        # merge prokka and counts
        merged = pd.merge(left=counts_df, right=prokka_df, how="left", left_on=COUNTS_HEADER[0], right_on=PROKKA_TSV_HEADER[1])

    logging.info("Parsing RefSeq file: %s" % refseq_tsv)
    refseq_df = get_valid_dataframe(refseq_tsv, REFSEQ_TSV_HEADER, sep="\t")

    if counts_tsv:
        # merge in refseq data
        merged = pd.merge(left=merged, right=refseq_df, how="left", left_on=COUNTS_HEADER[0], right_on=REFSEQ_TSV_HEADER[1])
    else:
        merged = pd.merge(left=prokka_df, right=refseq_df, how="left", left_on=PROKKA_TSV_HEADER[1], right_on=REFSEQ_TSV_HEADER[1])
    return merged


def merge_tables(prokka_tsv, refseq_tsv, output, counts_tsv=None):
    """

    Count data is a TSV formatted with a header:

        \b
        Geneid Chr Start  End Strand Length /path/example.bam
        orf1     1     1  500      +    500                50
        orf2     1   601  900      +    300               300
        orf3     1  1201 1500      +    300               200
    """
    df = do_merge(prokka_tsv, refseq_tsv, counts_tsv)
    if counts_tsv:
        df.to_csv(output, sep="\t", columns=MERGED_HEADER, index=False)
    else:
        df.to_csv(output, sep="\t", columns=ANNOTATE_HEADER, index=False)
