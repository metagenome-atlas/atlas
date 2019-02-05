#!/usr/bin/env python
import argparse
from snakemake.shell import shell
from snakemake.io import glob_wildcards



def predict_genes(genome,fasta,out_dir,log):

        fna = f"{out_dir}/{genome}.fna"
        faa = f"{out_dir}/{genome}.faa"
        gff = f"{out_dir}/{genome}.gff"

    shell(
        "prodigal -i {fasta} -o {gff} -d {fna} "
        "    -a {faa} -p meta -f gff 2> {log} "
        )


def predict_genes_genomes(in_dir,out_dir,log):

    path= os.path.join(in_dir,'{genome}.fasta')

    for genome in glob_wildcards(path).genome:
        predict_genes(genome,path.format(genome=genome),out_dir,log='log.txt')



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
