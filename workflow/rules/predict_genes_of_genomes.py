#!/usr/bin/env python3

import os, sys
import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logging.captureWarnings(True)


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



### New script

from Bio import SeqIO
import pyrodigal
from pathlib import Path


from multiprocessing import Pool
import itertools



def predict_genes(fasta_file, output_folder, translation_table=11, meta=True, name="infer") -> None:
    """Produces faa and fna file from fasta file using pyrodigal """

    if not name=="infer":
        # get name
        if fasta_file.endswith(".gz"): raise NotImplementedError("I have not implmented gziped fasta")
        name = Path(fasta_file).stem

    # define output files
    output_folder = Path(output_folder)
    faa_file= output_folder/ (name+".faa")
    fna_file= output_folder/ (name+".fna")


    orf_finder = pyrodigal.OrfFinder(meta=meta)

    with open(faa_file,"w") as faa_handle, open(fna_file,"w") as fna_handle:

        for contig in SeqIO.parse(fasta_file, "fasta"):

            genes = orf_finder.find_genes(bytes(contig.seq))


            genes.write_genes(fna_handle, sequence_id=contig.id)
            genes.write_translations(faa_handle, sequence_id=contig.id,translation_table=translation_table)


def predict_genes_genomes(input_dir, out_dir, threads):

    input_dir = Path(input_dir)
    out_dir = Path(out_dir)


    genomes_fastas= input_dir.glob("*[.fasta|.fn]*")

    out_dir.mkdir(exist_ok=True)
    
    pool = Pool(threads)
    pool.starmap(
        predict_genes,
        zip(genomes_fastas, itertools.repeat(out_dir),itertools.repeat(11)),
    )


if __name__ == "__main__":
    predict_genes_genomes(
        snakemake.input.dir,
        snakemake.output[0],
        int(snakemake.threads),
    )
