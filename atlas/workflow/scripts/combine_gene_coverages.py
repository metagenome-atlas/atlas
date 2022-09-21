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

import psutil

def measure_memory(write_log_entry=True):
    mem_uage = psutil.Process().memory_info().rss / (1024 * 1024)

    if write_log_entry:
        logging.info(f"The process is currelnty using {mem_uage: 7.0f} MB of RAM")

    return mem_uage


logging.info("Start")
measure_memory()

# read gene info

#gene_info= pd.read_table(input.info, index_col=0)
#gene_info.sort_index(inplace=True)
#gene_list= gene_info.index

#N_genes= gene_info.shape[0]

N_samples = len(snakemake.input.covstats)


output_dir = Path(snakemake.output[0])
output_dir.mkdir()


for i,cov_file in enumerate(snakemake.input.covstats):

    sample = os.path.split(cov_file)[-1].split("_")[0]
    logging.info(f"Read coverage file for sample {i+1}: {sample}")




    data = read_pileup_coverage(
        cov_file, coverage_measure="Avg_fold"
    )


    # transform index to int this should drastrically redruce memory
    data.index= pd.to_numeric(data.index.str[len("Gene"):].astype(int), downcast="integer")
    data.index.name = "GeneNr"
    # genes are not sorted
    data.sort_index(inplace=True)


    data["Avg_fold"] = pd.to_numeric(data.Avg_fold, downcast="float")
    data["Reads"] = pd.to_numeric(data.Reads, downcast="integer")


    output_file = output_dir /f"Sample={sample}"/"0.parq"
    output_file.parent.mkdir()

    data.to_parquet(output_file,index=True)

    # delete interminate data and release mem
    del data
    gc.collect()

    
    current_mem_uage= measure_memory()



logging.info("Finished")
