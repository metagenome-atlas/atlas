#!/usr/bin/env python3

# python 3.5 without f strings

import argparse
import os, shutil, sys
import uuid
import itertools
from glob import glob
from snakemake.shell import shell
from snakemake.io import glob_wildcards
from multiprocessing import Pool


def predict_genes(genome, fasta, out_dir, log):

    fna = "{}/{}.fna".format(out_dir, genome)
    faa = "{}/{}.faa".format(out_dir, genome)
    gff = "{}/{}.gff".format(out_dir, genome)

    shell('printf "{genome}:\n" > {log}'.format(genome=genome, log=log))
    shell(
        "prodigal -i {fasta} -o {gff} -d {fna} -a {faa} -p sinlge -c -m -f gff 2>> {log} ".format(
            fasta=fasta, log=log, gff=gff, fna=fna, faa=faa
        )
    )
    shell('printf "\n" >> {log}'.format(log=log))


def predict_genes_genomes(input_dir, out_dir, log, threads):
    genomes_fastas = glob(os.path.join(input_dir, "*.fasta"))

    os.makedirs(out_dir, exist_ok=True)

    temp_log_dir = os.path.join(os.path.dirname(log), "tmp_" + uuid.uuid4().hex)
    os.makedirs(temp_log_dir, exist_ok=False)

    genome_names = []
    log_names = []
    for fasta in genomes_fastas:
        genome_name = os.path.splitext(os.path.split(fasta)[-1])[0]
        genome_names.append(genome_name)
        log_names.append(os.path.join(temp_log_dir, genome_name + ".prodigal.tmp"))

    pool = Pool(threads)
    pool.starmap(
        predict_genes,
        zip(genome_names, genomes_fastas, itertools.repeat(out_dir), log_names),
    )

    # cat in python
    with open(log, "ab") as f_out:
        for logfile in log_names:
            with open(logfile, "rb") as f_in:
                shutil.copyfileobj(f_in, f_out)

    shell("rm -r {temp_log_dir}".format(temp_log_dir=temp_log_dir))


if __name__ == "__main__":
    try:
        log = open(snakemake.log[0], "w")
        sys.stderr = log
        sys.stdout = log

        predict_genes_genomes(
            snakemake.input.dir,
            snakemake.output[0],
            snakemake.log[0],
            int(snakemake.threads),
        )

    except NameError:
        p = argparse.ArgumentParser()
        p.add_argument("--input-dir", required=True)
        p.add_argument("--out-dir", required=True)
        p.add_argument("--log", required=True)
        p.add_argument("--threads", required=False, default=1, type=int)
        args = vars(p.parse_args())
        predict_genes_genomes(**args)
