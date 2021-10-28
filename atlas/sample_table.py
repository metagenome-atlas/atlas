import pandas as pd

ADDITIONAL_SAMPLEFILE_HEADERS = []  # ,'Contigs']

import logging

logger = logging.getLogger(__file__)


def validate_sample_table(sampleTable):

    Expected_Headers = ["BinGroup"] + ADDITIONAL_SAMPLEFILE_HEADERS
    for h in Expected_Headers:
        if not (h in sampleTable.columns):
            logger.error(f"expect '{h}' to be found in samples.tsv")
            exit(1)
        elif sampleTable[h].isnull().any():
            logger.error(f"Found empty values in the sample table column '{h}'")
            exit(1)

    if not sampleTable.index.is_unique:
        duplicated_samples = ", ".join(sampleTable.index.duplicated())
        logger.error(
            f"Expect Samples to be unique. Found {duplicated_samples} more than once"
        )
        exit(1)


    if sampleTable.index.str.match("^\d").any():
        logger.warning(
            f"Sample names shouldn't start with a digit. This can lead to incompatibilities.\n {list(sampleTable.index)}"
        )

    if sampleTable.index.str.contains('_').any():
        logger.warning(
            f"Sample names shouldn't contain underscores. This can lead to incompatibilities. \n {list(sampleTable.index)}"
        )
    if sampleTable.index.str.count('-').max()>1:
        logger.warning(
            f"Sample names shouldn't have more than one hypon '-'. This can lead to incompatibilities.\n {list(sampleTable.index)}"
        )


def load_sample_table(sample_table="samples.tsv"):

    sampleTable = pd.read_csv(sample_table, index_col=0, sep="\t")
    validate_sample_table(sampleTable)
    return sampleTable
