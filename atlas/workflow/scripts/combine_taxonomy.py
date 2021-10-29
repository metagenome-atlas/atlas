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

#### Begining of scripts

import pandas as pd
import numpy as np
from utils.taxonomy import tax2table

from glob import glob

gtdb_classify_folder = snakemake.input.folder

taxonomy_files = glob(f"{gtdb_classify_folder}/gtdbtk.*.summary.tsv")

N_taxonomy_files = len(taxonomy_files)
logging.info(f"Found {N_taxonomy_files} gtdb taxonomy files.")

if (0 == N_taxonomy_files) or (N_taxonomy_files > 2):

    raise Exception(
        f"Found {N_taxonomy_files} number of taxonomy files 'gtdbtk.*.summary.tsv' in {gtdb_classify_folder} expect 1 or 2."
    )


DT = pd.concat([pd.read_table(file, index_col=0) for file in taxonomy_files], axis=0)

DT.to_csv(snakemake.output.combined)

Tax = tax2table(DT.classification, remove_prefix=True)
Tax.to_csv(snakemake.output.taxonomy, sep="\t")
