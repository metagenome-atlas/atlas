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

    logger.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


# Install exception handler
sys.excepthook = handle_exception

#### Begining of scripts

import pandas as pd
import numpy as np
from utils.taxonomy import tax2table

from glob import glob
gtdb_classify_folder= snakemake.input.folder

taxonomy_files= glob(f"{gtdb_classify_folder}/gtdbtk.*.summary.tsv")

logging.INFO(f"Found {len(taxonomy_files)} gtdb taxonomy files.")
if 0==len(taxonomy_files) or len(taxonomy_files) <2:

    raise Exception(f"Found inadequate number of taxonomy files 'gtdbtk.*.summary.tsv' in {gtdb_classify_folder}"
                    )





DT= pd.concat([pd.read_table(file,index_col=0) for file in taxonomy_files],axis=0)

DT.to_csv(snakemake.output.combined)

Tax= tax2table(DT.classification, remove_prefix=True)
Tax.to_csv(snakemake.output.taxonomy,sep='\t')
