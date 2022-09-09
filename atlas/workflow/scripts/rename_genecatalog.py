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


map_genenr = pd.read_csv(snakemake.input.rep2genenr, index_col=0, sep="\t").squeeze()


# from gene Nr to gene name
rep2gene = geneNr_to_string(map_genenr)

logging.info(
    f"Collect and rename representative genes according to:\n {rep2gene.head()}"
)

assert rep2gene.shape[0] > 0

name_iterator = rep2gene.iteritems()
i=0
with open(snakemake.output[0], "w") as fout:
    with open(snakemake.input.fasta,"r") as fin:
        for line in fin:
            if line[0]==">":
                gene_name=line[1:].split(' ')[0]

                rep, gene_id = next(name_iterator)

                if ((i% 100)==0 ):
                    assert rep == gene_name, "genes are not in same order"

                fout.write(f">{gene_id} {gene_name}\n")
                i+=1
            else:
                fout.write(line)


   
