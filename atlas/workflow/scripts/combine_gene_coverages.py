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

# prepare snakemake.output tables
combined_cov = {}
combined_N_reads = {}


measure_memory()

N_samples = len(snakemake.input.covstats)
for i,cov_file in enumerate(snakemake.input.covstats):

    sample = os.path.split(cov_file)[-1].split("_")[0]





    data = read_pileup_coverage(
        cov_file, coverage_measure="Avg_fold"
    )

    # genes are not sorted
    # data.sort_index(inplace=True)


    combined_cov[sample] = pd.to_numeric(data.Avg_fold, downcast="float").astype(pd.SparseDtype(float, fill_value=0))
    combined_N_reads[sample] = pd.to_numeric(data.Reads, downcast="integer").astype(pd.SparseDtype(int, fill_value=0))

    # delete interminate data and release mem
    del data
    gc.collect()

    logging.info(f"Read coverage file for sample {i+1}: {sample}")
    current_mem_uage= measure_memory()
    estimated_max_mem = current_mem_uage/(i+1)*(N_samples+1)/1024

    logging.info(f"Estimated max mem is {estimated_max_mem:5.0f} GB")




# merge N reads
logging.info("Concatenate raw reads")
combined_N_reads = pd.concat(combined_N_reads,axis=1,sort=True,copy=False).fillna(0)
# give index nice name
combined_N_reads.index.name = "GeneNr"

measure_memory()



# Store as parquet
# add index so that it can be read in R

logging.info("Write first table")

from scipy.sparse import save_npz
save_npz("Genecatalog/counts/Nmapped_reads.npz", combined_N_reads.sparse.to_coo())
pd.Series(combined_N_reads.columns).to_csv("Genecatalog/counts/header.tsv",sep="\t", header=False,index=False)
pd.Series(combined_N_reads.index).to_csv("Genecatalog/counts/index.tsv",sep="\t", header=False,index=False)

#combined_N_reads.reset_index().to_hdf(,key='data',complevel=7,mode='w')
del combined_N_reads
gc.collect()
measure_memory()
# .sparse.to_coo() 

logging.info("Concatenate coverage data")
combined_cov = pd.concat(combined_cov,axis=1,sort=True,copy=False).fillna(0)
combined_cov.index.name = "GeneNr"

#combined_cov.reset_index().to_hdf(snakemake.output[0],key='data',complevel=7,mode='w')
save_npz("Genecatalog/counts/median_coverage.npz", combined_cov.sparse.to_coo())
del combined_cov
