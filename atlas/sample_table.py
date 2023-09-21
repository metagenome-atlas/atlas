import pandas as pd

import logging

logger = logging.getLogger(__file__)


def validate_sample_table(sampleTable):
    Expected_Headers = ["BinGroup"]
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
        logger.error(
            f"Sample names shouldn't start with a digit. This can lead to incompatibilities.\n {list(sampleTable.index)}"
        )
        exit(1)

    if sampleTable.index.str.contains("_").any():
        logger.error(
            f"Sample names shouldn't contain underscores. This can lead to incompatibilities. \n {list(sampleTable.index)}"
        )
        exit(1)

    if sampleTable.index.str.count("-").max() > 1:
        logger.error(
            f"Sample names shouldn't have more than one hypo '-'. This can lead to incompatibilities.\n {list(sampleTable.index)}"
        )
        exit(1)

    ### Validate BinGroup

    if sampleTable.BinGroup.isnull().any():
        logger.warning(f"Found empty values in the sample table column 'BinGroup'")

    if sampleTable.BinGroup.str.contains("_").any():
        logger.error(
            f"BinGroup names shouldn't contain underscores. This can lead to incompatibilities. \n {list(sampleTable.BinGroup)}"
        )
        exit(1)

    if sampleTable.BinGroup.str.contains("-").any():
        logger.error(
            f"BinGroup names shouldn't contain hypos '-'. This can lead to incompatibilities.\n {list(sampleTable.BinGroup)}"
        )
        exit(1)


def load_sample_table(sample_table="samples.tsv"):
    sampleTable = pd.read_csv(sample_table, index_col=0, sep="\t")
    validate_sample_table(sampleTable)
    return sampleTable


class BinGroupSizeError(Exception):
    """
    Exception with Bingroupsize
    """

    def __init__(self, message):
        super(BinGroupSizeError, self).__init__(message)


def validate_bingroup_size_cobinning(sampleTable, logger):
    """
    Validate that the bingroups are not too large, nor too small for co-binning.

    e.g. vamb and SemiBin
    """

    bin_group_sizes = sampleTable.BinGroup.value_counts()

    if bin_group_sizes.max() > 180:
        logger.warning(
            f"Found a bin group with more than 180 samples. This might lead to memory issues. \n {bin_group_sizes}"
        )

    if bin_group_sizes.min() < 10:
        logger.error(
            "If you want to use co-binning, you should have at least 5-10 samples per bin group. \n"
        )
        BinGroupSizeError("BinGroup too small")


def validate_bingroup_size_metabat(sampleTable, logger):
    bin_group_sizes = sampleTable.BinGroup.value_counts()

    max_bin_group_size = bin_group_sizes.max()

    warn_message = (
        "Co-binning with metabat uses cross-mapping which scales quadratically."
        f"You have a bingroup with {max_bin_group_size} samples, which already leads to {max_bin_group_size*max_bin_group_size} cross-mappings."
    )

    if max_bin_group_size > 50:
        logger.error(
            warn_message
            + "This is too much for metabat. Please use vamb, or SemiBin or split your samples into smaller groups."
        )
        BinGroupSizeError("BinGroup too large")

    if max_bin_group_size > 15:
        logger.warning(
            warn_message
            + "This might be too much for metabat. Consider using vamb, or SemiBin or split your samples into smaller groups."
        )

    elif max_bin_group_size == 1:
        logger.warning(
            "You have only one sample per bingroup. This doesn't use the co-abundance information."
        )


def validate_bingroup_size(sampleTable, config, logger):
    if config["final_binner"] == "DASTool":
        binners = config["binner"]

        logger.info(f"DASTool uses the folowing binners: {binners}")

        if ("vamb" in binners) or ("SemiBin" in binners):
            validate_bingroup_size_cobinning(sampleTable, logger)

        if "metabat" in binners:
            validate_bingroup_size_metabat(sampleTable, logger)

    elif config["final_binner"] == "metabat":
        validate_bingroup_size_metabat(sampleTable, logger)

    elif config["final_binner"] in ["vamb", "SemiBin"]:
        validate_bingroup_size_cobinning(sampleTable, logger)

    elif config["final_binner"] == "maxbin":
        logger.warning("maxbin Doesn't use coabundance for binning.")

    else:
        Exception(f"Unknown final binner: {config['final_binner']}")
