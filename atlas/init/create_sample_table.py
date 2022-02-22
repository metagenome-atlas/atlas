import os
import sys
import logging

logger = logging.getLogger(__file__)


import pandas as pd
import numpy as np
from collections import defaultdict


def splitext_ignore_gz(path):
    basename, ext = os.path.splitext(path)

    if ext == ".gz":
        basename, ext = os.path.splitext(basename)

    return basename, ext


def is_fastq_file(file, unzipped_extensions=[".fastq", ".fq"]):

    base_name, extension = splitext_ignore_gz(file)

    return extension in unzipped_extensions


def add_sample_to_table(sample_dict, sample_id, header, fastq):
    "Add fastq path to sample table, check if already in table"

    if (sample_id in sample_dict) and (header in sample_dict[sample_id]):

        logger.error(
            f"Duplicate sample {sample_id} was found with fraction {header};"
            f"\n Sample1: \n{sample_dict[sample_id]} \n"
            f"Sample2: {fastq}"
        )
        raise Exception()
    else:
        sample_dict[sample_id][header] = fastq


## global names
split_character = "infer"
is_paired = None


def infer_split_character(base_name):
    "Infer if fastq filename uses '_R1' '_1' to seperate filenames"

    global split_character, is_paired

    # infer split character if necessary only the first time.
    if (split_character is not None) and (split_character == "infer"):

        if ("_R1" in base_name) or ("_R2" in base_name):
            split_character = "_R"
            is_paired = True
        elif ("_1" in base_name) or ("_2" in base_name):
            split_character = "_"
            is_paired = True
        else:
            logger.warning(
                f"Could't find '_R1'/'_R2' or '_1'/'_2' in your filename {base_name}. Assume you have single-end reads."
            )
            split_character = None
            is_paired = False

        if split_character is not None:

            logger.info(
                f"I inferred that {split_character}1 and {split_character}2 distinguish paired end reads."
            )

    #    return split_character


def parse_file(full_path, sample_dict, sample_name=None):

    # only look at fastq files
    base_name, extension = splitext_ignore_gz(os.path.basename(full_path))

    # infer split character if necessary only the first time.
    infer_split_character(base_name)

    # pe reads
    if is_paired:

        if sample_name is None:
            sample_name = base_name.split(split_character)[0]

        if (split_character + "2") in base_name:
            add_sample_to_table(sample_dict, sample_name, "R2", full_path)
        elif (split_character + "1") in base_name:
            add_sample_to_table(sample_dict, sample_name, "R1", full_path)
        elif "_se" in base_name:
            logger.info("Found se reads. I ignore them.")
        else:

            logger.error(
                f"Did't find '{split_character}1' or  "
                f"'{split_character}2' in fastq {sample_name} : {full_path}"
                "Ignore file."
            )

    #  se reads
    else:

        if sample_name is None:
            sample_name = base_name

        add_sample_to_table(sample_dict, sample_name, "R1", full_path)


def parse_folder(base_folder, subfolder, sample_dict):

    files = os.listdir(os.path.join(base_folder, subfolder))

    fastq_files = [f for f in files if is_fastq_file(f)]

    N_files = len(fastq_files)

    if N_files > 3:

        # logger.info("You have more than three files per folder. Probably you want to merge them.")
        logger.critical(
            "You have more than three files per subfolder. I assume they come from different lanes. "
            "I cannot automatically merge them, yo you have to do it yourself. "
            " If you want tis feature. Write in https://github.com/metagenome-atlas/atlas/issues/389 "
        )
        raise NotImplementedError("Merging of files is not yet implemented")

    else:

        for file in fastq_files:
            parse_file(
                os.path.join(base_folder, subfolder, file),
                sample_dict,
                sample_name=subfolder,
            )


def get_samples_from_fastq(path, fraction_split_character=split_character):
    """
    creates table sampleID R1 R2 with the absolute paths of fastq files in a given folder
    """

    global split_character
    split_character = fraction_split_character

    sample_dict = defaultdict(dict)

    # list files and subfolder of fastq folder

    _, subfolders, files = next(os.walk(path))

    abs_path = os.path.abspath(path)

    # parse files
    for f in files:
        if is_fastq_file(f):
            parse_file(full_path=os.path.join(abs_path, f), sample_dict=sample_dict)

    # parse subfolder
    if len(subfolders) > 0:
        logger.info(
            f"Found {len(subfolders)} subfolders. Check if I find fastq files inside. Use the the subfolder as sample_names "
        )

        for subf in subfolders:
            parse_folder(abs_path, subf, sample_dict=sample_dict)

    # Create dataframe
    sample_df = pd.DataFrame(sample_dict).T.sort_index()

    if sample_df.isnull().any().any():
        logger.error(f"Missing files:\n\n {sample_df}")
        raise Exception()

    if sample_df.shape[0] == 0:
        logger.error(
            f"No files found in {path}\n"
            "    I'm looking for files with .fq or .fastq extension. "
        )
        raise Exception()

    logger.info(f"Found {sample_df.shape[0]} samples")

    return sample_df


def simplify_sample_names(sample_df):
    """Check if index of dataframe correspond to criteria and try to simplify them"""

    assert sample_df.index.is_unique

    sample_name_df = (
        sample_df.index.str.split("[_, ,-]", expand=True)
        .to_frame()
        .reset_index(drop=True)
    )

    # Index need simplification
    if sample_name_df.shape[1] > 1:
        logger.debug("Index need simplification")
        logger.debug(sample_name_df)
        sample_df["Full_Name"] = sample_df.index

        # First column is unique
        if sample_name_df[0].is_unique:
            logger.debug("First column is unique")
            sample_df.index = sample_name_df[0]

        # first two columns are unique
        elif not any(sample_name_df.duplicated(subset=[0, 1])):
            logger.debug("First two columns are unique: ")

            sample_df.index = sample_name_df.apply(
                lambda row: "{0}-{1}".format(*row), axis=1
            )

        # cannt find unique sample ids
        else:

            logger.warning(
                "Didn't found a way to simplify sample names. "
                "I imoute sample_names simply as S1...\n"
                "Modify the 'samples.tsv' created by atlas to suit this condition."
                "Sample names should consist only out of letters and numbers and start with a letter. "
            )

            sample_df.index = [f"S{i}" for i in range(1, 1 + sample_df.shape[0])]

        logger.info(
            "I simplified the sample names from {before} -> {after}\n".format(
                before=sample_df["Full_Name"].iloc[0], after=sample_df.index[0]
            )
        )

    # Start with digit
    if sample_df.index.str.match("^\d").any():

        logger.warning(f"Sample names start with a digit. I prepend 'S' for all.")
        sample_df.index = "S" + sample_df.index


### Testing


import os
import shutil

# create test folder structure


def create_test_fastq_files(
    output_folder,
    paired=True,
    fraction_indicator="_R",
    extension=".fastq.gz",
    subfolders=False,
    extra_file_name="L1_EXP5",
    samples=[f"sample{i}" for i in range(1, 4)],
):

    """Creates a folder (with subfolder) and touches fastq files for the test of atlas inti"""

    os.makedirs(output_folder)

    if paired:
        fractions = [f"{fraction_indicator}{r}" for r in range(1, 3)]
    else:
        fractions = [""]

    for s in samples:

        if subfolders:
            sample_split_character = "/"
        else:
            sample_split_character = "_"

        for fraction in fractions:

            fname = f"{output_folder}/{s}{sample_split_character}{extra_file_name}_{fraction}{extension}"
            os.makedirs(os.path.dirname(fname), exist_ok=True)

            open(fname, "w").close()


def test_table_creation(
    should_fail=False,
    samples=["sample1", "sample2", "sample3"],
    expected_samples=None,
    **kws,
):

    fastq_folder = "test_fastq_dir"
    if os.path.exists(fastq_folder):
        shutil.rmtree(fastq_folder)

    create_test_fastq_files(fastq_folder, samples=samples, **kws)

    try:
        sample_df = get_samples_from_fastq(fastq_folder)
        simplify_sample_names(sample_df)

    except Exception as e:

        if not should_fail:

            raise Exception("Test failed but shouldn't") from e

    if not should_fail:

        # check if samples are as expected
        if expected_samples is None:
            expected_samples = samples

        assert all(
            sample_df.index == expected_samples
        ), f"Samples not as expected {sample_df.index.values} != {expected_samples}"


def run_tests():

    for paired in [True, False]:

        test_table_creation(paired=paired)
        test_table_creation(
            samples=["S-ab-1", "S-ab-2"], expected_samples=["S1", "S2"], paired=paired
        )
        test_table_creation(samples=["S-1", "S-2"], paired=paired)
        test_table_creation(
            samples=["1", "2"], expected_samples=["S1", "S2"], paired=paired
        )

        test_table_creation(extension=".fq", paired=paired)
        test_table_creation(subfolders=True, paired=paired)

        test_table_creation(extra_file_name="_1F", paired=paired)
        test_table_creation(extension=".fasta", should_fail=True, paired=paired)

    test_table_creation(fraction_indicator="_")
    test_table_creation(extra_file_name="_1F", fraction_indicator="_", should_fail=True)


if __name__ == "__main__":
    run_tests()
