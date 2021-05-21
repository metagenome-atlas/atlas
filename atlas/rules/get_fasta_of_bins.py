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

    CA = pd.read_csv(cluster_attribution, header=None, index_col=1, sep="\t")

    assert CA.shape[1] == 1, "File should have only two columns " + cluster_attribution
    CA = CA.iloc[:, 0]
    CA.index = CA.index.astype("str")
    # exclude cluster 0 which is unclustered at least for metabat
    CA = CA.loc[CA != "0"]

    contigs = SeqIO.to_dict(SeqIO.parse(contigs, "fasta"))

    for id in CA.index.unique():
        bin_contig_names = CA.loc[id]
        out_file = os.path.join(out_folder, "{id}.fasta".format(id=id))
        if type(bin_contig_names) == str:
            warnings.warn("single contig bin Bin: " + out_file)
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

        with open(snakemake.log[0], "w") as log:
            sys.stderr = sys.stdout = log

            get_fasta_of_bins(
                snakemake.input.cluster_attribution,
                snakemake.input.contigs,
                snakemake.output[0],
            )
