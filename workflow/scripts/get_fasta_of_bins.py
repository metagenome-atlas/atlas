#! /usr/bin/env python


import sys, os
import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.DEBUG,
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


def get_fasta_of_bins(cluster_attribution, contigs_file, out_folder):
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

    contig_fasta_dict = SeqIO.index(contigs_file, "fasta")

    assert len(contig_fasta_dict) > 0, "No contigs in your fasta"

    unique_bins = CA.Bin.unique()

    assert len(unique_bins) >= 1, "No bins found"

    for binid in unique_bins:
        bin_contig_names = CA.loc[CA.Bin == binid, "Contig"].tolist()
        out_file = os.path.join(out_folder, f"{binid}.fasta")

        assert (
            len(bin_contig_names) >= 1
        ), f"No contigs found for bin {binid} in {cluster_attribution}"

        if len(bin_contig_names) == 1:
            warnings.warn(f"single contig bin Bin : {binid} {bin_contig_names}")

        logging.debug(f"Found {len(bin_contig_names)} contigs {bin_contig_names}")

        fasta_contigs = [contig_fasta_dict[c] for c in bin_contig_names]
        SeqIO.write(fasta_contigs, out_file, "fasta")


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
