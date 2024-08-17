# from ..color_logger import logger
import logging

logger = logging.getLogger(__file__)
import pandas as pd


Expected_library_values = {
    "library_selection": "RANDOM",
    "library_strategy": "WGS",
    "library_source": "METAGENOMIC",
    "platform": "ILLUMINA",
}

# Values that should be the same for all runs from the same sample
Expected_same_values = ["experiment_accession", "model", "library_name"]


SAMPLE_NAME_COLUMN = "sample_accession"


def load_and_validate_runinfo_table(path="RunInfo.csv"):
    RunTable = pd.read_csv(path, index_col=0)

    # validate sra table
    format_error = False

    # check if all headers are present
    Expected_headers = [
        "library_layout",
        "library_source",
        "library_selection",
        "library_strategy",
        "sample_accession",
    ]
    for header in Expected_headers:
        if not header in RunTable.columns:
            logger.error(f"Didn't found expected header {header}")
            format_error = True

    # set to string all expected headers
    RunTable[Expected_headers] = RunTable[Expected_headers].astype(str)

    if not all(RunTable.index.str[1:3] == "RR"):
        logger.error("Expect runs as index, e.g. [E,S,D]RR000")
        format_error = True

    assert (
        SAMPLE_NAME_COLUMN == "sample_accession"
    ), "This code is not tested for other sample names"
    if not (RunTable.sample_accession.str[1:3] == "RS").all():
        logger.error("sample_accession should start with [E,S,D]RS")
        format_error = True

    if not RunTable.library_layout.isin(["PAIRED", "SINGLE"]).all():
        logger.error("library_layout should be 'PAIRED' or 'SINGLE'")
        format_error = True

    if format_error:
        logger.error("RunTable {} is not valid. Abort.".format(path))
        exit(1)

    return RunTable


def filter_runinfo(RunTable, ignore_paired=False):
    logger.info(
        f"Start with {RunTable.shape[0]} runs from {RunTable[SAMPLE_NAME_COLUMN].unique().shape[0]} samples"
    )

    # Filter out reads that are not metagenomics

    for key in ["library_source"]:
        Nruns_before = RunTable.shape[0]
        All_values = RunTable[key].unique()
        RunTable = RunTable.loc[RunTable[key] == Expected_library_values[key]]

        Difference = Nruns_before - RunTable.shape[0]

        if Difference > 0:
            logger.info(
                f"Runs have the following values for {key}: {', '.join(All_values)}\n"
                f"Select only runs {key} == {Expected_library_values[key]}, "
                f"Filtered out {Difference} runs"
            )

    for key in ["library_selection", "library_strategy"]:
        Nruns_before = RunTable.shape[0]
        All_values = RunTable[key].unique()
        if any(RunTable[key] != Expected_library_values[key]):
            logger.warning(
                f"Runs have the following values for {key}: {', '.join(All_values)}\n"
                f"Usually I expect {key} == {Expected_library_values[key]} "
            )

    # Handle single end reads if mixed

    if ("PAIRED" in RunTable.library_layout) and ("SINGLE" in RunTable.library_layout):
        N_library_layout = RunTable.library_layout.value_counts()

        logger.info(
            f"Run table contains {N_library_layout['SINGLE']} single-end "
            f"and {N_library_layout['PAIRED']} paired-end libraries. "
        )

        if ignore_paired:
            logger.info(f"I drop {N_library_layout['PAIRED']} paired end libraries")
            RunTable = RunTable.query("library_layout == 'SINGLE'")

        else:
            logger.warning(f"I drop {N_library_layout['SINGLE']} single end libraries")

            RunTable = RunTable.query("library_layout == 'PAIRED'")

    # Illumina or not

    if not RunTable.platform.isin(["ILLUMINA"]).all():
        platforms = ", ".join(RunTable.platform.unique())

        logger.warning(
            f"Your samples are sequenced on the following platform: {platforms}\n"
            "I don't know how well Atlas handles non-illumina reads.\n"
            "If you have long-reads, specify them via a the longreads, column in the sample table."
        )

    # Final
    if RunTable.shape[0] > 0:
        logger.info(
            f"Selected {RunTable.shape[0]} runs from {RunTable[SAMPLE_NAME_COLUMN].unique().shape[0]} samples"
        )

    else:
        logger.critical("No runs left after filtering. Abort.")
        exit(1)

    return RunTable


def validate_merging_runinfo(path):
    RunTable = load_and_validate_runinfo_table(path)

    # If each run is from a different SAMPLE_NAME_COLUMN, merging is not necessary
    if RunTable.shape[0] == RunTable[SAMPLE_NAME_COLUMN].unique().shape[0]:
        return RunTable

    # Cannot merge if different platforms
    problematic_samples = []
    for sample, df in RunTable.groupby(SAMPLE_NAME_COLUMN):
        if not all(df.platform == df.platform.iloc[0]):
            problematic_samples.append(sample)

    if len(problematic_samples) > 0:
        logger.error(
            f"You attempt to merge runs from the same sample. "
            f"But for {len(problematic_samples)} samples the runs are sequenced with different platforms and shouldn't be merged.\n"
            f"Please resolve the ambiguity in the table {path} and rerun the command.\n"
        )

        exit(1)

    # Warn if samples are not identical values if expected the same

    for key in Expected_same_values:
        problematic_samples = []

        if key not in RunTable.columns:
            logger.warning(
                f"Didn't found column {key} in RunTable at {path}, don't check for identical values"
            )
        else:

            for sample, df in RunTable.groupby(SAMPLE_NAME_COLUMN):
                if not all(df[key] == df[key].iloc[0]):
                    problematic_samples.append(sample)

        if len(problematic_samples) > 0:
            if len(problematic_samples) > 5:
                problematic_samples_list = " ".join(problematic_samples[:3] + ["..."])
            else:
                problematic_samples_list = " ".join(problematic_samples)

                logger.warning(
                    "You attempt to merge runs from the same sample. "
                    f"But for {len(problematic_samples)} samples the runs have different {key}: {problematic_samples_list}\n"
                    f"You can modify the table {path} and rerun the command.\n"
                )

    logger.info("I will automatically merge runs from the same sample.")

    return RunTable


def get_run_ids_for_sample(run_table, sample):

    return run_table.query(f"{SAMPLE_NAME_COLUMN} == '{sample}'").index.tolist()


def get_all_sample_names(run_table):
    return run_table[SAMPLE_NAME_COLUMN].unique().tolist()
