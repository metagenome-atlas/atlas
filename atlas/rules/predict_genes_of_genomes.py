#!/usr/bin/env python3

# python 3.5 without f strings

import argparse
import os
from snakemake.shell import shell
from snakemake.io import glob_wildcards

def predict_genes(genome,fasta,out_dir,log):

    fna = "{}/{}.fna".format(out_dir,genome)
    faa = "{}/{}.faa".format(out_dir,genome)
    gff = "{}/{}.gff".format(out_dir,genome)

    shell("prodigal -i {fasta} -o {gff} -d {fna} -a {faa} -p meta -f gff 2>> {log} ".format(
        fasta=fasta, log=log,gff=gff,fna=fna,faa=faa)
          )



def predict_genes_genomes(in_dir,out_dir,log):

    path= os.path.join(in_dir,'{genome}.fasta')
    os.makedirs(out_dir,exist_ok=True)

    for genome in glob_wildcards(path).genome:
        predict_genes(genome,path.format(genome=genome),out_dir,log)



if __name__ == "__main__":
    try:

        predict_genes_genomes(
            snakemake.input[0],
            snakemake.output[0],
            snakemake.log[0]
        )

    except NameError:
        p = argparse.ArgumentParser()
        p.add_argument("--in-dir")
        p.add_argument("--out-dir")
        args = vars(p.parse_args())
        predict_genes_genomes(**args)
