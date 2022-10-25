import pandas as pd
import os
import numpy as np


def simplify_index(index):
    "assumes single index are path of files, removes extesnion and dirname"

    path = index[0]
    filename, extension = os.path.splitext(path)

    if extension == ".gz":
        extension = os.path.splitext(filename)[-1] + extension

    dir_name = os.path.dirname(path) + "/"

    return pd.Index(index).str.replace(extension, "").str.replace(dir_name, "")


def load_ani_table_(dist_file, header=None, simplify_names=False):

    F = pd.read_csv(dist_file, sep="\t", header=None, index_col=[0, 1])

    if header is not None:
        F.columns = header
    F.index.names = ["Genome1", "Genome2"]

    if simplify_names:

        F.index = pd.MultiIndex(
            levels=[
                simplify_index(F.index.levels[0]),
                simplify_index(F.index.levels[1]),
            ],
            codes=F.index.codes,
        )

    return F

def convert_percentages(df):
    """Convet all columns with strings and % at the end to percentages"""
    for col in df.columns:
        if df.dtypes[col] == "object":
            if df[col].iloc[0].endswith("%"):
                df.loc[:, col] = df[col].str.rstrip("%").astype("float") / 100.0
                
def load_bbsketch(dist_file, format=3, simplify_names=True):
    """ reads output of sendsketch.sh
        format=3 [query,ref,ANI..]
        format=2 Table for one query
            parses parameters in first line returns df,params
    """

    if format == 3:
        bbs = pd.read_csv(dist_file, index_col=[0, 1], sep="\t")
        bbs.index.names = ["Genome1", "Genome2"]
        if (bbs.QTaxID == -1).all():
            bbs.drop(["QTaxID", "RTaxID"], axis=1, inplace=True)

        bbs["Identity"] = bbs.iloc[:, 0] / 100.0

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

        f = open(dist_file)
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


def load_fastani(dist_file, simplify_names=True):
    """Loads fastANI output calculates overlap.
    Outputs a table with ['Genome1','Genome2','Identity','Nmapped','Ntotal','Overlap' ] in header"""

    F = load_ani_table_(
        dist_file, ["Identity", "Nmapped", "Ntotal"], simplify_names=simplify_names
    )
    F.loc[:, "Overlap"] = F.Nmapped.values / F.Ntotal.values
    F["Identity"] /= 100.0

    return F


def load_mash(dist_file, simplify_names=True):
    """Loads fastANI output calculates overlap.
    Outputs a table with ['Genome1','Genome2','Distance','Pvalue','Fraction','Identity'] in header"""
    F = load_ani_table_(
        dist_file, ["Distance", "Pvalue", "Fraction"], simplify_names=simplify_names
    )

    for col in ["Distance", "Pvalue"]:
        F[col] = pd.to_numeric(F[col], downcast="float")

    F["Identity"] = 1 - F.Distance
    # parse fraction
    F["Nmapped"] = pd.to_numeric( F.Fraction.map(lambda s: int(s.split("/")[0])), downcast = 'integer')
    F["Ntotal"] = pd.to_numeric( F.Fraction.map(lambda s: int(s.split("/")[1])), downcast = 'integer')
    F["Fraction"] = pd.to_numeric( F.Nmapped / F.Ntotal, downcast = 'float')

    return F


def load_parquet(parquet_file):

    M= pd.read_parquet(parquet_file,columns=["Distance"])
    M['Identity']= 1-M.Distance
    return M

def load_bindash(dist_file, simplify_names=True):
    """Loads bindash output.
    Outputs a table with
    ['Genome1','Genome2','Distance','Pvalue','Fraction','Nmapped','Ntotal','Identity']
    in header.

    Bindash tables are not necessarily simetrical.
    """
    F = load_ani_table_(
        dist_file, ["Distance", "Pvalue", "Fraction"], simplify_names=simplify_names
    )

    F["Nmapped"] = F.Fraction.map(lambda s: int(s.split("/")[0])).astype(int)
    F["Ntotal"] = F.Fraction.map(lambda s: int(s.split("/")[1])).astype(int)
    F["Fraction"] = F.Nmapped / F.Ntotal
    F["Identity"] = 1 - F.Distance

    return F


def load_mummer(dist_file):

    M = pd.read_csv(dist_file, sep="\t", index_col=[0, 1])
    M["Identity"] = M.ANI
    return M


def load_minimap(dist_file):

    M = pd.read_csv(dist_file, sep="\t", index_col=[0, 1])
    assert "Identity" in M.columns
    return M
