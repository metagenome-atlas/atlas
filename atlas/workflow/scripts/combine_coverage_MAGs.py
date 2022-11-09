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
import os, gc
from utils.parsers_bbmap import read_coverage_binned, combine_coverages


contig2genome = pd.read_csv(
    snakemake.input.contig2genome, header=None, index_col=0, sep="\t"
).iloc[:, 0]


# sum counts
logging.info("Loading counts and coverage per contig")

combined_cov, Counts_contigs = combine_coverages(
    snakemake.input.coverage_files, snakemake.params.samples
)

combined_cov = combined_cov.T

combined_cov.insert(
    0, "Genome", value=pd.Categorical(contig2genome.loc[combined_cov.index].values)
)

logging.info(f"Saving coverage to {snakemake.output.coverage_contigs}")

combined_cov.reset_index().to_parquet(snakemake.output.coverage_contigs)

logging.info("Sum counts per genome")

Counts_genome = Counts_contigs.groupby(contig2genome, axis=1).sum().T
Counts_genome.index.name = "Sample"

logging.info(f"Saving counts to {snakemake.output.counts}")

Counts_genome.reset_index().to_parquet(snakemake.output.counts)
del Counts_genome, combined_cov, Counts_contigs
gc.collect()

# Binned coverage
logging.info("Loading binned coverage")
binCov = {}
for i, cov_file in enumerate(snakemake.input.binned_coverage_files):

    sample = snakemake.params.samples[i]

    binCov[sample] = read_coverage_binned(cov_file)

binCov = pd.DataFrame.from_dict(binCov)

logging.info("Add genome information to it")
binCov.insert(
    0, "Genome", value=pd.Categorical(contig2genome.loc[binCov.index.get_level_values(0)].values)
)

gc.collect()
logging.info(f"Saving combined binCov to {snakemake.output.binned_cov}")
binCov.reset_index().to_parquet(snakemake.output.binned_cov)

# Median coverage
logging.info("Calculate median coverage")
Median_abund = binCov.groupby("Genome").median().T
del binCov
gc.collect()
logging.info(f"Saving mediuan coverage {snakemake.output.median_abund}")
Median_abund.reset_index().to_parquet(snakemake.output.median_abund)
