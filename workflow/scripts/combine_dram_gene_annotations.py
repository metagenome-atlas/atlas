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


from pathlib import Path
import numpy as np
import pandas as pd
from collections import defaultdict

db_columns = {
    "kegg": ["ko_id", "kegg_hit"],
    "peptidase": [
        "peptidase_id",
        "peptidase_family",
        "peptidase_hit",
        "peptidase_RBH",
        "peptidase_identity",
        "peptidase_bitScore",
        "peptidase_eVal",
    ],
    "pfam": ["pfam_hits"],
    "cazy": ["cazy_ids", "cazy_hits", "cazy_subfam_ec", "cazy_best_hit"],
    # "heme": ["heme_regulatory_motif_count"],
}

Tables = defaultdict(list)

for file in snakemake.input:
    df = pd.read_csv(file, index_col=0, sep="\t")

    # drop un-annotated genes
    df = df.query("rank!='E'")

    # change index from 'subset1_Gene111' ->  simply 'Gene111'
    # Gene name to nr
    df.index = (
        df.index.str.split("_", n=1, expand=True)
        .get_level_values(1)
        .str[len("Gene") :]
        .astype(np.int64)
    )
    df.index.name = "GeneNr"

    # select columns, drop na rows and append to list
    for db in db_columns:
        cols = db_columns[db]

        if not df.columns.intersection(cols).empty:
            Tables[db].append(df[cols].dropna(axis=0, how="all"))

    del df

out_dir = Path(snakemake.output[0])
out_dir.mkdir()

for db in Tables:
    combined = pd.concat(Tables[db], axis=0)

    combined.sort_index(inplace=True)

    combined.reset_index().to_parquet(out_dir / (db + ".parquet"))
