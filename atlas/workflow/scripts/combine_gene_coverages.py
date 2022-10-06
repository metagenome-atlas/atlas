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

import psutil


def measure_memory(write_log_entry=True):
    mem_uage = psutil.Process().memory_info().rss / (1024 * 1024)

    if write_log_entry:
        logging.info(f"The process is currelnty using {mem_uage: 7.0f} MB of RAM")

    return mem_uage


logging.info("Start")
measure_memory()

# read gene info

# gene_info= pd.read_table(input.info, index_col=0)
# gene_info.sort_index(inplace=True)
# gene_list= gene_info.index

# N_genes= gene_info.shape[0]

N_samples = len(snakemake.input.covstats)

# prepare snakemake.output tables
combined_cov = {}
combined_N_reads = {}


for i, sample in enumerate(snakemake.params.samples):

    cov_file = snakemake.input.covstats[i]

    data = read_pileup_coverage(cov_file, coverage_measure="Median_fold")

    # transform index to int this should drastrically redruce memory
    data.index = data.index.str[len("Gene") :].astype(int)

    # genes are not sorted
    # data.sort_index(inplace=True)

    combined_cov[sample] = pd.to_numeric(data.Median_fold, downcast="float")
    combined_N_reads[sample] = pd.to_numeric(data.Reads, downcast="integer")

    # delete interminate data and release mem
    del data
    gc.collect()

    logging.info(f"Read coverage file for sample {i+1}: {sample}")
    current_mem_uage = measure_memory()
    estimated_max_mem = current_mem_uage / (i + 1) * (N_samples + 1) / 1024

    logging.info(f"Estimated max mem is {estimated_max_mem:5.0f} GB")


# merge N reads
logging.info("Concatenate raw reads")
combined_N_reads = pd.concat(combined_N_reads, axis=1, sort=True, copy=False).fillna(0)
# give index nice name
combined_N_reads.index.name = "GeneNr"

measure_memory()


# Store as parquet
# add index so that it can be read in R

logging.info("Write first table")

combined_N_reads.reset_index().to_parquet(snakemake.output[1])
del combined_N_reads
gc.collect()
measure_memory()


logging.info("Concatenate coverage data")
combined_cov = pd.concat(combined_cov, axis=1, sort=True, copy=False).fillna(0)
combined_cov.index.name = "GeneNr"

combined_cov.reset_index().to_parquet(snakemake.output[0])

del combined_cov
