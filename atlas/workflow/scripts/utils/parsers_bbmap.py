import pandas as pd
import os


def parse_comments(file, comment="#", sep="\t", expect_one_value=True):
    "parse comments at begin of file #Avg: 123"
    Parsed = {}
    with open(file) as f:
        line = f.readline()

        while line.startswith(comment):
            line_values = line[1:].strip().split(sep)
            name = line_values[0]
            if name[-1] == ":":
                name = name[:-1]
            values = line_values[1:]

            if len(values) == 1:
                Parsed[name] = values[0]
            elif not expect_one_value:
                Parsed[name] = values

            line = f.readline()

    if len(Parsed) == 0:
        raise Exception(
            "Couldn't parse values from file {} with comment char {} and sep '{}' ".format(
                file, comment, sep
            )
        )
    return Parsed


def read_coverage_binned(covarage_binned_file):

    binned_coverage = pd.read_csv(
        covarage_binned_file,
        sep="\t",
        skiprows=2,
        index_col=[0, 2],
        usecols=[0, 1, 2],
    )
    binned_coverage.index.names = ["Contig", "Position"]
    binned_coverage.columns = ["Coverage"]

    return binned_coverage.Coverage


def read_pileup_coverage(pileup_file, coverage_measure="Median_fold", other_columns=[]):
    """
    Parse output of pileup.sh
    #ID             Avg_fold        Length  Ref_GC  Covered_percent Covered_bases   Plus_reads      Minus_reads     Read_GC Median_fold     Std_Dev
    sample1_0       126.9594        80077   0.2537  100.0000        80077           46301           46301           0.2569  128             29.59
    sample1_1       124.7803        63666   0.2644  100.0000        63666           36186           36186           0.2678  128             20.08
    sample1_2       125.6420        62268   0.2576  100.0000        62268           35642           35642           0.2614  128             19.89
    sample1_3       126.0232        54480   0.2664  100.0000        54480           31283           31283           0.2686  128             17.25

    """

    data = pd.read_csv(
        pileup_file,
        index_col=0,
        sep="\t",
        usecols=["#ID", coverage_measure, "Plus_reads", "Minus_reads"]
        + list(other_columns),
    )
    data.index.names = ["Contig"]
    data.loc[data[coverage_measure] < 0, coverage_measure] = 0
    data.eval("Reads = Plus_reads + Minus_reads", inplace=True)
    data.drop(["Plus_reads", "Minus_reads"], axis=1, inplace=True)

    return data


def combine_coverages(coverage_files, sample_names, coverage_measure="Median_fold"):
    """
    Combines the coverage files from different samples
    Args:
        coverage_files: bunch of coverage_files produced with pileup.sh from the bbmap package
        sample_names: sample names associated with the coverage_files

    Output:
        combined_cov:   pandas dataframe of samples x contigs for coverage
        combined_N_reads:   pandas dataframe of samples x contigs for number of mapped reads
    """

    combined_cov = {}
    combined_N_reads = {}

    assert len(coverage_files) == len(sample_names)

    for i in range(len(coverage_files)):

        sample = sample_names[i]

        data = read_pileup_coverage(
            coverage_files[i], coverage_measure=coverage_measure
        )

        combined_cov[sample] = data[coverage_measure]
        combined_N_reads[sample] = data.Reads

    combined_cov = pd.DataFrame(combined_cov).T
    combined_N_reads = pd.DataFrame(combined_N_reads).T

    return combined_cov, combined_N_reads


def parse_bbmap_log_file(log_file):
    """
    parses a bbmap log file (paired, single end or both (bbrwap))
    returns number of used and mapped reads.
    This is the sum or se + R1 + R2 reads
    """
    N_results = 0

    mapped = 0
    used = 0

    with open(log_file) as f:
        for line in f:
            if line.startswith("mapped:"):
                # mapped:                  85.4123%        27804399        85.4133%         2774727540
                try:
                    mapped_reads = line.strip().split(":")[1].split()[1]
                except ValueError as e:
                    raise Exception(f"Error parsing line:\n{line}\n") from e
                mapped += int(mapped_reads)

            elif line.startswith("Reads Used:"):
                # Reads Used:             65106274        (6496839447 bases)

                try:
                    used_reads = line.strip().split(":")[1].split()[0]
                except ValueError as e:
                    raise Exception(f"Error parsing line:\n{line}\n") from e

                used += int(used_reads)

            elif "Results" in line:
                N_results += 1
                assert N_results <= 2, "you have more than one library in this log file"

        if used == 0:
            raise IOError(
                f"I couldn't parse the log file, probably the job didn't finish properly: {log_file}"
            )

        assert (
            used >= mapped
        ), "something is wrong, you have more than 100% mapped reads?"

        return used, mapped


# def parse_simple_log_file(log_file, keyword, expect_one_value=True):
#     content = open(log_file).read()
#     pos = content.find(keyword)
#     if pos == -1:
#         raise Exception("Didn't find {} in file:\n\n{}".format(keyword, log_file))
#
#     else:
#         if expect_one_value:
#             return content[pos:].split()[2]
#
#         else:
#             return content[pos:].split()[2:]
