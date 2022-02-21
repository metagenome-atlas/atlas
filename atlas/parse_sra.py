from ..color_logger import logger
import pandas as pd


def load_and_validate_sra_table(path):



    Expected_library_values= {#"LibrarySelection": "RANDOM",
                              "LibraryStrategy": "WGS" ,
                              "LibrarySource": "METAGENOMIC"
                             }

    for key in Expected_library_values:
        if not all( RunTable[key] == Expected_library_values[key] ):
            logger.critical(f"'{key}' should be '{Expected_library_values[key]}' for all reads")

            exit(1)


    


    Nlanes= pd.crosstab(RunTable.BioSample,RunTable.LibraryLayout)
    assert ~ Nlanes.isnull().any().any()
    assert (Nlanes.diff(axis=1).iloc[:,1]==0).all(), f"Not all samples have equal lanes {Nlanes}"

    return RunTable


def get_sra_fractions():

    LibraryLayouts= list( RunTable.LibraryLayout.unique())

    if 'PAIRED' in LibraryLayouts:
        Fractions = ['1','2']
    else:
        Fractions = []
    if 'SINGLE' in LibraryLayouts:
        Fractions += ['se']

    assert len(Fractions)>0

    return Fractions



def validate_runinfo_table(RunTable):

    # validate sra table
    Expected_headers= ['LibraryLayout','LibrarySource','LibrarySelection','LibraryStrategy','BioSample']
    for header in Expected_headers:
        if not header in RunTable.columns:

            logger.error(f"Didn't found expected header {header} in sra table {path}"
                            )
            exit(1)


    assert all(RunTable.index.str[1:2]=="R"), "Expect runs as index"

    assert RunTable.BioSample.str.startswith('SAM').all(), "BioSample should start with 'SAM'"

    assert RunTable.LibraryLayout.isin(['PAIRED','SINGLE']).all() , "LibraryLayout should be paired or single"



def filter_runinfo(RunInfo):


    Expected_library_values= { "LibrarySelection": "RANDOM",
                               "LibraryStrategy": "WGS" ,
                               "LibrarySource": "METAGENOMIC"
    }


                             

    # Filter out reads that are not metagenomics

    for key in Expected_library_values:

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

        logger.info("I drop single end libraries")

        RunTable = RunTable.query("LibraryLayout == 'PAIRED'")


    # Final 

    logger.info(f"Selected {RunTable.shape[0]} runs from {RunTable.BioSample.unique().shape[0]} samples")


    return RunTable




    