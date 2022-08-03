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
from utils.parsers_bbmap import read_bbsplit_bincov
import gc

# Binned coverage

binCov = {}
for i, cov_file in enumerate(snakemake.input.binned_coverage_files):

    sample = snakemake.params.samples[i]
    logging.info(f"Reading bincov_file for sample {sample}")

    binCov[sample] = read_bbsplit_bincov(cov_file)

binCov = pd.DataFrame.from_dict(binCov)

gc.collect()
logging.info(f"Saving combined binCov to {snakemake.output.binCov}")
binCov.reset_index().to_parquet(snakemake.output.binned_cov)


# Median coverage
logging.info("Calculate median coverage")
Median_abund = binCov.groupby(level=0).median().T
del binCov
Median_abund.reset_index().to_parquet(snakemake.output.median_abund)


del Median_abund
gc.collect()
