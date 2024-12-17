import os, sys
import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logging.captureWarnings(True)


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

#### Beginning of scripts

import pandas as pd
from utils.parsers import read_checkm_output


def main(samples, completeness_files, taxonomy_files, bin_table):
    sample_data = {}
    div = {}

    df = pd.DataFrame()

    for i, sample in enumerate(samples):
        sample_data = read_checkm_output(
            taxonomy_table=taxonomy_files[i], completness_table=completeness_files[i]
        )
        sample_data["Sample"] = sample

        df = df.append(sample_data)

    df.to_csv(bin_table, sep="\t")


if __name__ == "__main__":
    main(
        samples=snakemake.params.samples,
        taxonomy_files=snakemake.input.taxonomy_files,
        completeness_files=snakemake.input.completeness_files,
        bin_table=snakemake.output.bin_table,
    )
