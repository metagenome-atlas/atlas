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

from utils.fasta import parse_fasta_headers
from utils.utils import gen_names_for_range
from glob import glob
import pandas as pd


fasta_files = glob(f"{snakemake.input[0]}/*{snakemake.params.extension}")

if len(fasta_files) > 0:

    Bin_names = gen_names_for_range(
        N=len(fasta_files), prefix=f"{snakemake.wildcards.sample}_SemiBin_"
    )

    mappings = []

    for bin_name, fasta in zip(Bin_names, fasta_files):
        contigs = parse_fasta_headers(fasta)

        mappings.append(pd.Series(data=bin_name, index=contigs))

    pd.concat(mappings, axis=0).to_csv(
        snakemake.output[0], sep="\t", header=False, index=True
    )

else:

    logging.warning(
        f"No bins found in {snakemake.input[0]} add longest contig as bin to make atlas continue."
    )

    with open(snakemake.output[0], "w") as outf:
        outf.write("{sample}_0\t{sample}_SemiBin_1\n".format(**snakemake.wildcards))
