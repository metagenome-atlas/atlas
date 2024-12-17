import pandas as pd
import networkx as nx
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc
import numpy as np

import logging

logger = logging.getLogger(__name__)


def verify_expected_range(value, min_range, max_range, name="value"):
    if (value < min_range) or (value > max_range):
        raise Exception(
            f"{name} should be in range [{min_range},{max_range}] but is {value}"
        )


def load_bbsketch(dist_file, format=3, simplify_names=True):
    """reads output of sendsketch.sh
    format=3 [query,ref,ANI..]
    format=2 Table for one query
        parses parameters in first line returns df,params
    """

    if format == 3:
        bbs = pd.read_csv(dist_file, index_col=[0, 1], sep="\t")
        bbs.index.names = ["Genome1", "Genome2"]
        if (bbs.QTaxID == -1).all():
            bbs.drop(["QTaxID", "RTaxID"], axis=1, inplace=True)

        bbs["ANI"] = bbs.iloc[:, 0] / 100.0

        if "SSU" in bbs:
            bbs["SSU"] = bbs.SSU.replace(".", np.nan)

        if simplify_names:
            bbs.index = pd.MultiIndex(
                levels=[
                    simplify_index(bbs.index.levels[0]),
                    simplify_index(bbs.index.levels[1]),
                ],
                codes=bbs.index.codes,
            )

        return bbs
    elif format == 2:
        f = open(send_sketch_file)
        f.readline()  # trash empty line
        comment_line = f.readline().strip()
        params = dict(key_value.split(":") for key_value in comment_line.split("\t"))

        df = pd.read_csv(f, sep="\t")

        convert_percentages(df)

        return df, params
    else:
        raise NotImplementedError(
            "I don't know how to parse other formats than 2,3 of bbsketch"
        )


def define_aligned_fraction(skani_df):
    "Define aligned fraction by max"

    skani_df["Align_fraction"] = skani_df[
        ["Align_fraction_ref", "Align_fraction_query"]
    ].max(axis=1)


def load_skani(parquet_file):
    M = pd.read_parquet(
        parquet_file,
        columns=["Ref", "Query", "ANI", "Align_fraction_ref", "Align_fraction_query"],
    )

    # make fraction not percentages
    for col in ["ANI", "Align_fraction_ref", "Align_fraction_query"]:
        M.eval(f" {col} = {col} / 100", inplace=True)

    define_aligned_fraction(M)

    assert "Align_fraction" in M.columns, M.columns

    M.set_index(["Ref", "Query"], inplace=True)

    return M


def load_bindash(dist_file, simplify_names=True):
    """Loads bindash output.
    Outputs a table with
    ['Genome1','Genome2','Distance','Pvalue','Fraction','Nmapped','Ntotal','ANI']
    in header.

    Bindash tables are not necessarily symmetrical.
    """
    F = load_ani_table_(
        dist_file, ["Distance", "Pvalue", "Fraction"], simplify_names=simplify_names
    )

    F["Nmapped"] = F.Fraction.map(lambda s: int(s.split("/")[0])).astype(int)
    F["Ntotal"] = F.Fraction.map(lambda s: int(s.split("/")[1])).astype(int)
    F["Fraction"] = F.Nmapped / F.Ntotal
    F["ANI"] = 1 - F.Distance

    return F


def to_graph(F, attributes=None, **kws):
    df = F.copy()

    df["Genome1"] = df.index.get_level_values(0)
    df["Genome2"] = df.index.get_level_values(1)

    G = nx.from_pandas_edgelist(df, "Genome1", "Genome2", attributes, **kws)

    return G


def best_genome_from_table(Grouping, quality_score):
    Mapping = pd.Series(index=Grouping.index)

    for group in Grouping.unique():
        genomes = Grouping.index[Grouping == group]
        representative = quality_score.loc[genomes].idxmax()
        Mapping.loc[genomes] = representative

    return Mapping


def clustermap(DistanceMatrix, linkage_method="average", **kws):
    import seaborn as sns
    import scipy.spatial as sp, scipy.cluster.hierarchy as hc

    linkage = hc.linkage(sp.distance.squareform(DistanceMatrix), method=linkage_method)

    cg = sns.clustermap(
        1 - DistanceMatrix, row_linkage=linkage, col_linkage=linkage, **kws
    )

    return cg


def pairewise2matrix(M, fillna=np.nan):
    """
    This functions turns a pairewise genome ANI table [genome1, genome2, column...]
    In to a matrix [genome1 genome2] of the values of column.
    When ANI values are symmetrical (with minimal error),
    usually only one halve of NxN possibilities values are calculated.

    During the process missing values are inputted with 0

    Diagonal values are set to 1

    """
    assert type(M) == pd.Series

    if (M < 0).any():
        raise Exception("Some Id values are < 0")

    if M.isnull().any():
        raise Exception("Some Id values are NA")

    # check if not zeros
    n_zeros = (M == 0).sum()
    if (n_zeros > 0) and (fillna != 0):
        logger.warning(
            f"{n_zeros} of id values are zero, they will be replaced by {fillna}"
        )

    ID = M.unstack()

    all_indexes = ID.index.union(ID.columns)
    ID = ID.reindex(index=all_indexes, columns=all_indexes)
    ID = ID.fillna(0)
    ID = ID + ID.T
    ID.values[np.eye(ID.shape[0], dtype=bool)] = 1

    n_zeros = (ID == 0).sum().sum()
    if n_zeros > 0:
        logger.debug(
            f"Impute {n_zeros} ({n_zeros/ (ID.shape[0] * ID.shape[1])*100:.2}%) values with {fillna}"
        )

    if fillna == 0:
        return ID
    else:
        return ID.replace(0, fillna)


def group_species_linkage(M, threshold=0.95, fillna=None, linkage_method="average"):
    """
    M is a series of ANI
    """
    assert type(M) == pd.Series
    verify_expected_range(threshold, 0.3, 1, "clustering threshold")

    verify_expected_range(M.max(), 0.05, 1, "ANI max")
    verify_expected_range(M.min(), 0.05, 1, "ANI min")

    cutoff = 1 - threshold

    if fillna is None:
        fillna = M.min() * 0.65

    Dist = 1 - pairewise2matrix(M, fillna=fillna)

    linkage = hc.linkage(sp.distance.squareform(Dist), method=linkage_method)
    labels = pd.Series(
        hc.fcluster(linkage, cutoff, criterion="distance"), index=Dist.index
    )

    return labels
