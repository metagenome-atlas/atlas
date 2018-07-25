#! /bin/env python
import pandas as pd
from Bio import SeqIO

import argparse, os, shutils



def get_fasta_of_bins(cluster_attribution,contigs,out_prefix):
    """
        Creates individual fasta files for each bin using the contigs fasta and the cluster attribution.

        input:
        - cluster attribution file:   tab seperated file of "contig_fasta_header    bin"
        - contigs:                    fasta file of contigs
        - out_prefix:                 output_prefix for bin fastas  {output_prefix}.{binid}.fasta

    """
    # create outdir
    out_folder= os.path.dirname(out_prefix)
    if os.path.exists(out_folder): shutils.rmtree(out_folder)
        os.makedirs(out_folder)


    CA = pd.read_table(cluster_attribution,header=None)

    contigs = SeqIO.to_dict(SeqIO.parse(contigs,'fasta'))

    for id in CA.index.unique():

        bin_contigs = [contigs[c] for c in CA.loc[id]]

        out_file = out_prefix+ srt(id) + '.fasta'

        SeqIO.write(bin_contigs,out_file,'fasta')


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--cluster-attribution")
    p.add_argument("--contigs")
    p.add_argument("--out-prefix")
    args = vars(p.parse_args())

    get_fasta_of_bins(**args)
