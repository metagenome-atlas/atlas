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


# start of script
import argparse
import os, sys
import shutil
import warnings

import pandas as pd
from Bio import SeqIO


def get_fasta_of_bins(cluster_attribution, contigs, out_folder):
    """
    Creates individual fasta files for each bin using the contigs fasta and the cluster attribution.

    input:
    - cluster attribution file:   tab seperated file of "contig_fasta_header    bin"
    - contigs:                    fasta file of contigs
    - out_prefix:                 output_prefix for bin fastas  {out_folder}/{binid}.fasta
    """
    # create outdir
    if os.path.exists(out_folder):
        shutil.rmtree(out_folder)
    os.makedirs(out_folder)

    CA = pd.read_csv(cluster_attribution, header=None, sep="\t", dtype=str)

    assert CA.shape[1] == 2, "File should have only two columns " + cluster_attribution

    CA.columns = ["Contig", "Bin"]

    # assert that Contig is unique
    assert CA.Contig.is_unique, (
        f"First column of file {cluster_attribution} should be contigs, hence unique"
        f"I got\n{CA.head()}"
    )

    contigs = SeqIO.to_dict(SeqIO.parse(contigs, "fasta"))

    for binid in CA.Bin.unique():
        bin_contig_names = CA.query('Bin == "@binid"').index.tolist()
        out_file = os.path.join(out_folder, f"{binid}.fasta")

        if type(bin_contig_names) == str:
            warnings.warn(f"single contig bin Bin : {binid} {bin_contig_names}")
            bin_contig_names = [bin_contig_names]

        bin_contigs = [contigs[c] for c in bin_contig_names]
        SeqIO.write(bin_contigs, out_file, "fasta")


if __name__ == "__main__":
    if "snakemake" not in globals():
        p = argparse.ArgumentParser()
        p.add_argument("--cluster-attribution")
        p.add_argument("--contigs")
        p.add_argument("--out-folder")
        args = vars(p.parse_args())
        get_fasta_of_bins(**args)
    else:
        get_fasta_of_bins(
            snakemake.input.cluster_attribution,
            snakemake.input.contigs,
            snakemake.output[0],
        )
