#! /usr/bin/env python3


import sys, os
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

annotation_file = snakemake.input.annotations
module_output_table = snakemake.output[0]

from mag_annotator.database_handler import DatabaseHandler
from mag_annotator.summarize_genomes import build_module_net, make_module_coverage_frame

annotations = pd.read_csv(annotation_file, sep="\t", index_col=0)


# get db_locs and read in dbs
database_handler = DatabaseHandler(logger=logging, config_loc=snakemake.input.config)


if "module_step_form" not in database_handler.config["dram_sheets"]:
    raise ValueError(
        "Module step form location must be set in order to summarize genomes"
    )

module_steps_form = pd.read_csv(
    database_handler.config["dram_sheets"]["module_step_form"], sep="\t"
)

all_module_nets = {
    module: build_module_net(module_df)
    for module, module_df in module_steps_form.groupby("module")
}

module_coverage_frame = make_module_coverage_frame(
    annotations, all_module_nets, groupby_column="fasta"
)

module_coverage_frame.to_csv(module_output_table, sep="\t")
