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

#### Begining of scripts

import pandas as pd
from utils.parsers import read_busco_output


def main(samples, completeness_files, bin_table):
    sample_data = {}
    div = {}

    df = pd.DataFrame()

    for i, sample in enumerate(samples):
        sample_data = read_busco_output(completeness_files[i])
        sample_data["Sample"] = sample

        df = df.append(sample_data)

    # remove missing

    failed_genomes = df.index[df.Dataset.str.lower().str.contains("run failed")]

    if len(failed_genomes) > 0:
        logging.warn(
            "Following genomes didn't pass BUSCO. I ignore them, because "
            "I think theas means they are too bad to be quantified:\n"
            f"{failed_genomes}"
        )

        df.loc[failed_genomes, ["Completeness", "Contamination", "Quality_score"]] = 0

    df.to_csv(bin_table, sep="\t")


if __name__ == "__main__":
    main(
        samples=snakemake.params.samples,
        completeness_files=snakemake.input.completeness_files,
        bin_table=snakemake.output.bin_table,
    )
