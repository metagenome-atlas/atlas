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


import pandas as pd
from utils.gene_scripts import geneNr_to_string


# Start
import pyfastx

map_genenr = pd.read_csv(snakemake.input.rep2genenr, index_col=0, sep="\t").squeeze()


# from gene Nr to gene name
rep2gene = geneNr_to_string(map_genenr)

logging.info(
    f"Collect and rename representative genes according to:\n {rep2gene.head()}"
)

assert rep2gene.shape[0] > 0

fasta_parser = pyfastx.Fasta(snakemake.input.fasta, build_index=True)


with open(snakemake.output[0], "w") as fout:
    for orf, gene_name in rep2gene.iteritems():

        fout.write(f">{gene_name} {orf}\n{fasta_parser[orf].seq}\n")

