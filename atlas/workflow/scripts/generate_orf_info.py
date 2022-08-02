#! /usr/bin/env python


import sys, os
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


## Start

import pandas as pd
import numpy as np

from utils import gene_scripts

# CLuterID    GeneID    empty third column
orf2gene = pd.read_csv(
    snakemake.input.cluster_attribution, header=None, sep="\t", usecols=[0, 1]
)

orf2gene.columns = ["ORF", "Representative"]

# split orf names in sample, contig_nr, and orf_nr
orf_info = gene_scripts.split_orf_to_index(orf2gene.ORF)

# rename representative

representative_names = orf2gene.Representative.unique()

map_names = pd.Series(
    index=representative_names,
    data=np.arange(1, len(representative_names) + 1, dtype=np.uint),
)


orf_info["GeneNr"] = orf2gene.Representative.map(map_names)


orf_info.to_parquet(snakemake.output.cluster_attribution)


# Save name of representatives
map_names.index.name = "Representative"
map_names.name = "GeneNr"
map_names.to_csv(snakemake.output.rep2genenr, sep="\t")
