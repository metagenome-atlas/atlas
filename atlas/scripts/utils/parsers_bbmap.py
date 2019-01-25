

import pandas as pd
import os

def read_coverage_binned(covarage_binned_file):
    return pd.read_table(covarage_binned_file,
                     skiprows=2,
                     index_col=[0,2],
                     usecols=[0,1,2],
                     squeeze=True)

def combine_coverages(coverage_files,sample_names,coverage_measure='Median_fold'):
    """
        Combines the coverage files from different samples
        Args:
            coverage_files: bunch of coverage_files produced with pileup.sh from the bbmap package
            sample_names: sample names associated with the coverage_files

        Output:
            combined_cov:   pandas dataframe of samples x contigs for coverage
            combined_N_reads:   pandas dataframe of samples x contigs for number of mapped reads
    """

    combined_cov={}
    combined_N_reads={}

    assert len(coverage_files)==len(sample_names)

    for i in range(len(coverage_files)):

        sample= sample_names[i]
        data= pd.read_table(coverage_files[i],index_col=0)

        data.loc[data[coverage_measure]<0,coverage_measure]=0
        combined_cov[sample]= data[coverage_measure]
        combined_N_reads[sample] = data.Plus_reads+data.Minus_reads

    combined_cov= pd.DataFrame(combined_cov).T
    combined_N_reads= pd.DataFrame(combined_N_reads).T

    return combined_cov, combined_N_reads
