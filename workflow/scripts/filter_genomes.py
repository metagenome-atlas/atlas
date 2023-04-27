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


import pandas as pd
from glob import glob
from numpy import log

from utils.parsers import load_quality


Q = load_quality(snakemake.input.quality)

stats = pd.read_csv(snakemake.input.stats, index_col=0, sep="\t")
stats["logN50"] = log(stats.N50)

# merge table but only shared Bins and non overlapping columns
Q = Q.join(stats.loc[Q.index, stats.columns.difference(Q.columns)])
del stats

n_all_bins = Q.shape[0]

filter_criteria = snakemake.params["filter_criteria"]
logging.info(f"Filter genomes according to criteria:\n {filter_criteria}")


Q = Q.query(filter_criteria)

logging.info(f"Retain {Q.shape[0]} genomes from {n_all_bins}")


## GUNC

if hasattr(snakemake.input, "gunc"):

    gunc = pd.read_table(snakemake.input.gunc, index_col=0)
    gunc = gunc.loc[Q.index]

    bad_genomes = gunc.index[gunc["pass.GUNC"] == False]
    logging.info(f"{len(bad_genomes)} Don't pass gunc filtering")

    Q.drop(bad_genomes, inplace=True)
else:
    logging.info(" Don't filter based on gunc")


if Q.shape[0] == 0:
    logging.error(
        f"No bins passed filtering criteria! Bad luck!. You might want to tweek the filtering criteria. Also check the {snakemake.input.quality}"
    )
    exit(1)

# output Q together with quality
Q.to_csv(snakemake.output.info, sep="\t")


# output quality for derepliation
D = Q.copy()

D.index.name = "genome"

D.columns = D.columns.str.lower()
# fasta extension is needed even if otherwise stated https://github.com/MrOlm/drep/issues/169
D.index += ".fasta"
D = D[["completeness", "contamination"]]
D.to_csv(snakemake.output.quality_for_derep)

# filter path genomes

F = pd.read_table(snakemake.input.paths, index_col=0).squeeze()

F = F.loc[Q.index].iloc[:, 0]
F.to_csv(snakemake.output.paths, index=False, header=False)
