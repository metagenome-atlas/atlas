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

import h5py

import psutil
def measure_memory(write_log_entry=True):
    mem_uage = psutil.Process().memory_info().rss / (1024 * 1024)

    if write_log_entry:
        logging.info(f"The process is currelnty using {mem_uage: 7.0f} MB of RAM")

    return mem_uage


logging.info("Start")
measure_memory()

N_samples = len(snakemake.input.covstats)

logging.info("Read gene info")

gene_info= pd.read_table(input.info, index_col=0)
gene_info.sort_index(inplace=True)
N_genes= gene_info.shape[0]
#gene_list= gene_info.index

del gene_info
gc.collect()

measure_memory()



logging.info("Open hdf files for writing")

gene_matrix_shape = (N_samples,N_genes)

with h5py.File(output.cov, 'w') as hdf_cov_file, h5py.File(output.counts, 'w') as hdf_counts_file: 


    combined_cov = hdf_cov_file.create_dataset('data', shape= gene_matrix_shape ,fillvalue=0, compression="gzip")
    combined_counts = hdf_cov_file.create_dataset('data', shape= gene_matrix_shape ,fillvalue=0, compression="gzip")

    # add Smaple names attribute
    sample_names = np.array(list(snakemake.params.samples)).astype("S")
    combined_cov.attrs['sample_names'] = sample_names
    combined_counts.attrs['sample_names'] = sample_names
    

    Summary= {}


    logging.info("Start reading files")
    initial_mem_uage = measure_memory()


    for i,sample in enumerate(snakemake.params.samples):

        sample_cov_file = snakemake.input.covstats[i]

        data = read_pileup_coverage(sample_cov_file, coverage_measure="Median_fold")

        # transform index to int this should drastrically redruce memory
        # data.index = data.index.str[len("Gene") :].astype(int)

        # I hope genes are sorted
        assert data.index.is_monotonic_increasing

        # downcast data 
        Median_fold = pd.to_numeric(data.Median_fold, downcast="float")
        Reads= pd.to_numeric(data.Reads, downcast="integer")
        

        # delete interminate data and release mem
        del data
        gc.collect()


        # get summary statistics
        logging.info("Extract Summary statistics")
        non_zero_coverage= Median_fold.loc[Median_fold>0]
        
        Summary[sample] = {"Sum_coverage" : Median_fold.sum(), 
                           "Total_counts" : Reads.sum(), 
                           "Non_zero_counts"   : (Reads>0).sum(),
                           "Non_zero_coverage"   : non_zero_coverage.shape[0],
                           "Max_coverage" : Median_fold.max(),
                           "Average_nonzero_coverage" : non_zero_coverage.mean(),
                           "Q1_nonzero_coverage" : non_zero_coverage.percentile(0.25),
                           "Median_nonzero_coverage" : non_zero_coverage.percentile(0.5),
                           "Q3_nonzero_coverage" : non_zero_coverage.percentile(0.75),
                           }
        del non_zero_coverage

        logging.info(Summary[sample])

        combined_cov[i,:] = Median_fold.values
        combined_counts[i,:] = Reads.values

    
        del Median_fold,Reads
        gc.collect()


        logging.info(f"Read coverage file for sample {i+1} / {N_samples}")
        current_mem_uage = measure_memory()
        estimated_max_mem = (current_mem_uage- initial_mem_uage) /(i + 1) * (N_samples + 1) 

        logging.info(f"Estimated max mem is {estimated_max_mem:5.0f} MB")


logging.info("All samples processed")
gc.collect()

logging.info("Save Summary")
pd.DataFrame(Summary).to_csv(snakemake.output.summary,sep='\t')