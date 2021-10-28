import pandas as pd
import numpy as np
import warnings

TAXONMIC_LEVELS = ["Domain", "phylum", "class", "order", "family", "genus", "species"]


def tax2table(Taxonomy_Series, split_character=";", remove_prefix=False):
    """
    Transforms (green_genes) taxonomy to a table
    Expect the following input format:
    d__Bacteria;p__Bacteroidota;c__Bacteroidia;f__
    Replaces empty values and can remove prefix 'c__'
    """

    # drop missing values
    if Taxonomy_Series.isnull().any():
        warnings.warn(
            "Some samples have no taxonomy asigned. Samples:\n"
            + ", ".join(Taxonomy_Series.index[Taxonomy_Series.isnull()].astype(str))
        )

    Tax = Taxonomy_Series.dropna().astype(str).str.split(split_character, expand=True)
    # Add headers as long as we have columns
    Tax.columns = TAXONMIC_LEVELS[: len(Tax.columns)]

    if remove_prefix:
        Tax = Tax.applymap(lambda s: s[3:], na_action="ignore").replace("", np.nan)
    else:
        Tax[Tax.applymap(len, na_action="ignore") == 3] = np.nan

    # add missing values again

    Tax = Tax.reindex(Taxonomy_Series.index)

    return Tax


def load_checkm_tax(taxonomy_file, remove_prefix=False):

    D = pd.read_table(taxonomy_file, index_col=0)

    checkmTax = tax2table(D["Taxonomy (contained)"], remove_prefix=remove_prefix)

    return checkmTax


def load_gtdb_tax(taxonomy_file, remove_prefix=False):

    D = pd.read_table(taxonomy_file, index_col=0)

    Tax = tax2table(D["classification"], remove_prefix=remove_prefix)

    return Tax
