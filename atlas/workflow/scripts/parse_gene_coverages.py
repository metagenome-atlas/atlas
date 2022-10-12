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
from pathlib import Path




logging.info("Start")


# read gene info

#gene_info= pd.read_table(input.info, index_col=0)
#gene_info.sort_index(inplace=True)
#gene_list= gene_info.index

#N_genes= gene_info.shape[0]




output_file = Path(snakemake.output[0])

sample = snakemake.wildcards.sample

logging.info(f"Read coverage file for sample {sample}")


data = read_pileup_coverage(
    snakemake.input[0], coverage_measure="Median_fold"
)


# transform index to int this should drastrically redruce memory
data.index= pd.to_numeric(data.index.str[len("Gene"):].astype(int), downcast="integer")
data.index.name = "GeneNr"
# genes are not sorted
data.sort_index(inplace=True)

# Downcasting creates inconsistent datatypes
#data["Avg_fold"] = pd.to_numeric(data.Avg_fold, downcast="float")
#data["Reads"] = pd.to_numeric(data.Reads, downcast="integer")


data.to_parquet(output_file,index=True)



logging.info("Finished")