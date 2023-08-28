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


from Bio import SeqIO

# Open the snakemake.output FASTA file and mapping table file for writing
with open(snakemake.output.fasta, "w") as output_handle, open(
    snakemake.output.mapping_table, "w"
) as mapping_table_handle:
    i = 1

    for record in SeqIO.parse(snakemake.input[0], "fasta"):
        if len(record) < snakemake.params.minlength:
            break

        old_name = record.id
        new_name = f"{snakemake.wildcards.sample}_{i}"
        record.id = new_name
        record.description = ""

        SeqIO.write(record, output_handle, "fasta")

        mapping_table_handle.write(f"{new_name}\t{old_name}\n")

        i += 1
