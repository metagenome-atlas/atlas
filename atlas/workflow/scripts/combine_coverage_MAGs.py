import os, sys
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
import os
from utils.parsers_bbmap import read_coverage_binned, combine_coverages


contig2genome = pd.read_csv(
    snakemake.input.contig2genome, header=None, index_col=0, sep="\t"
).iloc[:, 0]


# sum counts

combined_cov, Counts_contigs = combine_coverages(
    snakemake.input.covstats, snakemake.params.samples
)


Counts_genome = Counts_contigs.groupby(contig2genome, axis=1).sum().T
Counts_genome.index.name = "Sample"
Counts_genome.to_csv(snakemake.output.counts, sep="\t")


# Binned coverage

binCov = {}
for i, cov_file in enumerate(snakemake.input.binned_coverage_files):

    sample = snakemake.params.samples[i]

    binCov[sample] = read_coverage_binned(cov_file)

binCov = pd.DataFrame.from_dict(binCov)
binCov.to_csv(snakemake.output.binned_cov, sep="\t", compression="gzip")


# Median coverage

Median_abund = (
    binCov.groupby(contig2genome.loc[binCov.index.get_level_values(0)].values)
    .median()
    .T
)

Median_abund.to_csv(snakemake.output.median_abund, sep="\t")
