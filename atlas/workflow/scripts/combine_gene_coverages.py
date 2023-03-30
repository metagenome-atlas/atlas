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
import numpy as np
import pandas as pd
import gc, os


import h5py

import h5py

import psutil
def measure_memory(write_log_entry=True):
    mem_uage = psutil.Process().memory_info().rss / (1024 * 1024)

    if write_log_entry:
        logging.info(f"The process is currently using {mem_uage: 7.0f} MB of RAM")

    return mem_uage


logging.info("Start")
measure_memory()

N_samples = len(snakemake.input.covstats)

logging.info("Read gene info")

gene_info= pd.read_table(snakemake.input.info)

# Gene name is only first part of first column
gene_info.index=gene_info["#Name"].str.split(" ",n=1,expand=True)[0]
gene_info.index.name= "GeneName"
gene_info.drop("#Name",axis=1,inplace=True)

gene_info= pd.read_table(snakemake.input.info, index_col=0)
gene_info.sort_index(inplace=True)
N_genes= gene_info.shape[0]
#gene_list= gene_info.index

# Sort 
gene_info.sort_index(inplace=True)
N_genes= gene_info.shape[0]

gene_info[["Samples_nz_coverage","Samples_nz_counts","Sum_coverage","Max_coverage"]]=0

 
#gene_list= gene_info.index



logging.info("Open hdf files for writing")

gene_matrix_shape = (N_samples,N_genes)

with h5py.File(snakemake.output.cov, 'w') as hdf_cov_file, h5py.File(snakemake.output.counts, 'w') as hdf_counts_file: 


    combined_cov = hdf_cov_file.create_dataset('data', shape= gene_matrix_shape ,fillvalue=0, compression="gzip")
    combined_counts = hdf_counts_file.create_dataset('data', shape= gene_matrix_shape ,fillvalue=0, compression="gzip")

    # add Smaple names attribute
    sample_names = np.array(list(snakemake.params.samples)).astype("S")
    combined_cov.attrs['sample_names'] = sample_names
    combined_counts.attrs['sample_names'] = sample_names
    
    gc.collect()
    
    Summary= {}


    logging.info("Start reading files")
    initial_mem_uage = measure_memory()


    for i,sample in enumerate(snakemake.params.samples):


        logging.info(f"Read coverage file for sample {i+1} / {N_samples}")
        sample_cov_file = snakemake.input.covstats[i]

        

        data = pd.read_parquet(sample_cov_file, columns=["GeneName","Reads","Median_fold"]).set_index("GeneName")

        assert data.shape[0] == N_genes, f"I only have {data.shape[0]} /{N_genes} in the file {sample_cov_file}"

        # genes are not sorted :-()
        assert data.index.is_monotonic_increasing,f"data is not sorted by index in {sample_cov_file}"

        

        # downcast data 
        # median is int
        Median_fold = pd.to_numeric(data.Median_fold, downcast="integer")
        Reads= pd.to_numeric(data.Reads, downcast="integer")
        

        # delete interminate data and release mem
        del data



        # get summary statistics per sample
        logging.debug("Extract Summary statistics")
        
        Summary[sample] = {"Sum_coverage" : Median_fold.sum(), 
                           "Total_counts" : Reads.sum(), 
                           "Genes_nz_counts"   : (Reads>0).sum(),
                           "Genes_nz_coverage"   : (Median_fold>0).sum()
                           }

        # get gene wise stats
        gene_info["Samples_nz_counts"] += (Reads>0)*1 
        gene_info["Samples_nz_coverage"] += (Median_fold>0)*1 
        gene_info["Sum_coverage"] += Median_fold

        gene_info["Max_coverage"] = np.fmax(gene_info["Max_coverage"],Median_fold)

        combined_cov[i,:] = Median_fold.values
        combined_counts[i,:] = Reads.values

    
        del Median_fold,Reads
        gc.collect()


        

        current_mem_uage = measure_memory()






logging.info("All samples processed")
gc.collect()

logging.info("Save sample Summary")
pd.DataFrame(Summary).to_csv(snakemake.output.sample_info,sep='\t')


logging.info("Save gene Summary")

# downcast 
for col in gene_info.columns:

    if col == "GC":
        gene_info[col] = pd.to_numeric(gene_info[col], downcast="float")
    else:
        gene_info[col] = pd.to_numeric(gene_info[col], downcast="integer")

gene_info.reset_index().to_parquet(snakemake.output.gene_info)
