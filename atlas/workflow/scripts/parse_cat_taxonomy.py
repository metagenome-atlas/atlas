#!/usr/bin/env python3


import pandas as pd
import numpy as np


def main(taxonomy_file, tax_output):
    """
    parses output of CAT v4.3.3 removes duplicated classification and makes nice table of only taxonomy names
    """

    DT = pd.read_csv(taxonomy_file, sep="\t")

    DT = DT.fillna("NoStdName")
    # sort by number of classified levels
    DT["Nunclassified"] = (DT.loc[:, "superkingdom":] == "not classified").sum(1)
    DT = DT.sort_values(["# bin", "Nunclassified"])

    if DT.duplicated("# bin", keep=False).shape[0] > 0:
        print(
            "Found duplicated taxonomic annotation, drop duplicates, keep first with more levels of classification."
        )
        print(DT.loc[DT.duplicated("# bin", keep=False)].T)

    DT = DT.drop_duplicates(["# bin"], keep="first")
    DT.index = DT["# bin"]
    DT.drop(["Nunclassified"], axis=1, inplace=True)

    Tax = DT.loc[:, "superkingdom":]
    Tax = Tax.applymap(lambda s: s.split(":")[0])
    Tax_nan = Tax.replace("not classified", np.nan)

    Tax_nan.to_csv(tax_output, sep="\t")


if __name__ == "__main__":

    try:
        main(snakemake.input[0], snakemake.output[0])
    except NameError:
        import argparse

        p = argparse.ArgumentParser()
        p.add_argument("-i", "--input-taxonomy", dest="taxonomy_file")
        p.add_argument("-o", "--parsed-taxonomy", dest="tax_output")
        args = vars(p.parse_args())
        main(**args)
        rename_genomes(**args)
