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
from utils import gene_scripts

# if MAGs are renamed I need to obtain the old contig names
# otherwise not
if snakemake.params.renamed_contigs:
    contigs2bins = pd.read_csv(
        snakemake.input.contigs2bins, index_col=0, squeeze=False, sep="\t", header=None
    )

    contigs2bins.columns = ["Bin"]
    old2newID = pd.read_csv(
        snakemake.input.old2newID, index_col=0, squeeze=True, sep="\t"
    )

    contigs2genome = contigs2bins.join(old2newID, on="Bin").dropna().drop("Bin", axis=1)
else:
    contigs2genome = pd.read_csv(
        snakemake.input.contigs2mags, index_col=0, squeeze=False, sep="\t", header=None
    )
    contigs2genome.columns = ["MAG"]

# load orf_info
orf_info = pd.read_parquet(snakemake.input.orf_info)


# recreate Contig name `Sample_ContigNr` and Gene names `Gene0004`
orf_info["Contig"] = orf_info.Sample + "_" + orf_info.ContigNr.astype(str)
orf_info["Gene"] = gene_scripts.geneNr_to_string(orf_info.GeneNr)

# Join genomes on contig
orf_info = orf_info.join(contigs2genome, on="Contig")

# remove genes not on genomes
orf_info = orf_info.dropna(axis=0)


# count genes per genome in a matrix
gene2genome = pd.to_numeric(
    orf_info.groupby(["Gene", "MAG"]).size(), downcast="unsigned"
).unstack(fill_value=0)

# save as parquet
gene2genome.reset_index().to_parquet(snakemake.output[0])
