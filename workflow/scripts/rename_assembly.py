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

import gzip as gz
from Bio import SeqIO

# Open the snakemake.output FASTA file and mapping table file for writing
with open(snakemake.output.mapping_table, "w") as mapping_table_handle, gz.open(snakemake.output.fasta_gz, "wt") as out_fasta_gz, open(snakemake.output.fasta, "w") as out_fasta_uncompressed, gz.open(snakemake.input[0], "rt") as in_fasta :
    i = 1

    for record in SeqIO.parse(in_fasta, "fasta"):

        old_name = record.id
        new_name = f"{snakemake.wildcards.sample}_{i}"
        record.id = new_name
        record.description = ""

        SeqIO.write(record, out_fasta_gz, "fasta")
        SeqIO.write(record, out_fasta_uncompressed, "fasta")

        mapping_table_handle.write(f"{new_name}\t{old_name}\n")

        i += 1
    
    assert i > 1, "No sequences found in input FASTA file"
