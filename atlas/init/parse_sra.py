# from ..color_logger import logger
import logging
logger = logging.getLogger(__file__)
import pandas as pd



Expected_library_values= { "LibrarySelection": "RANDOM",
                            "LibraryStrategy": "WGS" ,
                            "LibrarySource": "METAGENOMIC",
                            "Platform": "ILLUMINA"
}


def load_and_validate_runinfo_table(path):

    RunTable = pd.read_csv(path, sep="\t", index_col=0)

    # validate sra table
    format_error= False

    # check if all headers are present
    Expected_headers= ['LibraryLayout','LibrarySource','LibrarySelection','LibraryStrategy','BioSample']
    for header in Expected_headers:
        if not header in RunTable.columns:

            logger.error(f"Didn't found expected header {header}")
            format_error= True


    if not all(RunTable.index.str[1:2]=="R"):
        logger.error("Expect runs as index to start with 'R'")
        format_error= True

    if not RunTable.BioSample.str.startswith('SAM').all():
        logger.error("BioSample should start with 'SAM'")
        format_error= True

    if not RunTable.LibraryLayout.isin(['PAIRED','SINGLE']).all():
        logger.error("LibraryLayout should be 'PAIRED' or 'SINGLE'")
        format_error= True

    if format_error:
        logger.error("RunTable {} is not valid. Abort.".format(path))
        exit(1)

    return RunTable



def filter_runinfo(RunTable, ignore_paired=False):


    # Filter out reads that are not metagenomics

    for key in ["LibrarySource","LibrarySelection"]:

        Nruns_before= RunTable.shape[0]
        RunTable = RunTable.loc[ RunTable[key] == Expected_library_values[key] ]
        
        Difference = Nruns_before - RunTable.shape[0]

        if Difference > 0:

            logger.info(f"Select only runs {key} == {Expected_library_values[key]},"
                        f" Filtered out {Difference} runs"
                        )

    


    # Handle single end reads if mixed

    if ('PAIRED' in RunTable.LibraryLayout) and ('SINGLE' in RunTable.LibraryLayout):


        N_library_layout= RunTable.LibraryLayout.value_counts()

        logger.info(f"Run table contains {N_library_layout['SINGLE']} single-end "
                     f"and {N_library_layout['PAIRED']} paired-end libraries. "
                     )

        if ignore_paired:
            logger.info(f"I drop {N_library_layout['PAIRED']} paired end libraries")
            RunTable = RunTable.query("LibraryLayout == 'SINGLE'")

        else:

            logger.warn(f"I drop {N_library_layout['SINGLE']} single end libraries")

            RunTable = RunTable.query("LibraryLayout == 'PAIRED'")


    # Illumina or not

    if not RunTable.Platform.isin(['ILLUMINA']).all():
        Platforms= ", ".join(RunTable.Platform.unique())
        
        logger.warn(f"Your samples are sequenced on the folowing platform: {Platforms}\n"
                    "I don't know how well Atlas handles non-illumina reads.\n"
                    "If you have long-reads, specify them via a the longreads, column in the sample table."
                    )


    # Final 

    logger.info(f"Selected {RunTable.shape[0]} runs from {RunTable.BioSample.unique().shape[0]} samples")


    return RunTable


def validate_merging_runinfo(path):

    RunTable = load_and_validate_runinfo_table(path)

    # If each run is from a different biosample, merging is not necessary
    if RunTable.shape[0] == RunTable.BioSample.unique().shape[0]:
        return RunTable
    

    

    # Cannot merge if different platforms
    problematic_samples=[]
    for sample, df in RunTable.groupby('BioSample'):
        if not all(df.Platform == df.Platform.iloc[0]):
            problematic_samples.append(sample)


    if len(problematic_samples) > 0:
        logger.error(f"You attemt to merge runs from the same sample. "
                     f"But for {len(problematic_samples)} samples the runs are sequenced with different platforms and should't be merged.\n" 
                     f"Please resolve the the abiguity in the table {path} and rerun the command.\n"
                    )

        exit(1)

    # Warn if samples are not identical for the follwing columns
    Expected_same_values = ["Experiment","Model","LibraryName"]
    for key in Expected_same_values:

        problematic_samples=[]
        for sample, df in RunTable.groupby('BioSample'):
            if not all(df[key] == df[key].iloc[0]):
                problematic_samples.append(sample)

        if len(problematic_samples) > 0:
            if len(problematic_samples)>5:
                problematic_samples_list = " ".join(problematic_samples[:3]+["..."])
            else:
                problematic_samples_list= " ".join(problematic_samples)

                logger.warn("You attemt to merge runs from the same sample. "
                            f"But for {len(problematic_samples)} samples the runs have different {key}: {problematic_samples_list}\n"
                            f"You can modify the table {path} and rerun the command.\n"
                            )

    logger.info("I will automatically merge runs from the same biosample.")

    return RunTable




    
        
                    

    