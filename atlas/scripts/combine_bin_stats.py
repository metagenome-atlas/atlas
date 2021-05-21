import os, sys

log = open(snakemake.log[0], "w")
sys.stderr = log
sys.stdout = log

import pandas as pd
from snakemake.utils import report


atlas_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.join(atlas_dir, "scripts"))
from utils.parsers_checkm import read_checkm_output


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
