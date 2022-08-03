#!/usr/bin/env python
import os, sys
import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logging.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


# Install exception handler
sys.excepthook = handle_exception

#### Begining of script

import pandas as pd
import gc, os
from utils.parsers_bbmap import read_pileup_coverage

# prepare snakemake.output tables
combined_cov = pd.DataFrame(
    dtype=int,
)
combined_N_reads = pd.DataFrame(
    dtype=float,
)


for cov_file in snakemake.input.covstats:

    sample = os.path.split(cov_file)[-1].split("_")[0]

    # print(f"read file for sample {sample}")

    if cov_file == snakemake.input.covstats[0]:
        other_columns = ["Length"]
    else:
        other_columns = []

    data = read_pileup_coverage(
        cov_file, coverage_measure="Avg_fold", other_columns=other_columns
    )

    # genes are not sorted
    data.sort_index(inplace=True)

    # add gene length to dataframe of counts
    if cov_file == snakemake.input.covstats[0]:
        combined_N_reads["Length"] = pd.to_numeric(
            data.Length, downcast="unsigned"
        )

    combined_cov[sample] = pd.to_numeric(data.Avg_fold, downcast="float")
    combined_N_reads[sample] = pd.to_numeric(data.Reads, downcast="unsigned")

    # delete interminate data and releace mem
    del data
    gc.collect()


# give index nice name
combined_cov.index.name = "GeneNr"
combined_N_reads.index.name = "GeneNr"

# Store as parquet
# add index so that it can be read in R
combined_N_reads.reset_index().to_parquet(snakemake.output[1])
del combined_N_reads
gc.collect()

combined_cov.reset_index().to_parquet(snakemake.output[0])
del combined_cov