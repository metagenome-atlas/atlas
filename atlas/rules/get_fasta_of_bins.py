import argparse
import os
import shutil
import warnings

import pandas as pd
from Bio import SeqIO


def get_fasta_of_bins(cluster_attribution, contigs, out_prefix):
    """
    Creates individual fasta files for each bin using the contigs fasta and the cluster attribution.

    input:
    - cluster attribution file:   tab seperated file of "contig_fasta_header    bin"
    - contigs:                    fasta file of contigs
    - out_prefix:                 output_prefix for bin fastas  {output_prefix}.{binid}.fasta
    """
    # create outdir
    out_folder = os.path.dirname(out_prefix)
    if os.path.exists(out_folder):
        shutil.rmtree(out_folder)
    os.makedirs(out_folder)

    CA = pd.read_table(cluster_attribution, header=None, index_col=1)

    assert CA.shape[1] == 1, "File should have only two columns " + cluster_attribution
    CA = CA.iloc[:, 0]
    CA.index = CA.index.astype("str")

    contigs = SeqIO.to_dict(SeqIO.parse(contigs, "fasta"))

    for id in CA.index.unique():
        bin_contig_names = CA.loc[id]
        out_file = "{prefix}.{id}.fasta".format(prefix=out_prefix, id=id)
        if type(bin_contig_names) == str:
            warnings.warn("single contig bin Bin: " + out_file)
            bin_contig_names = [bin_contig_names]
        bin_contigs = [contigs[c] for c in bin_contig_names]
        SeqIO.write(bin_contigs, out_file, "fasta")


if __name__ == "__main__":
    if snakemake is not None:
        get_fasta_of_bins(
            snakemake.input.cluster_attribution,
            snakemake.input.contigs,
            snakemake.params.prefix,
        )
    else:
        p = argparse.ArgumentParser()
        p.add_argument("--cluster-attribution")
        p.add_argument("--contigs")
        p.add_argument("--out-prefix")
        args = vars(p.parse_args())
        get_fasta_of_bins(**args)
