import os, sys, stat

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

    logger.error(
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
from snakemake.utils import report


def main(input_files, report_out):

    input_files = "\n".join(input_files)

    report_str = """


=============================================================
ATLAS_ - Dummy Report
=============================================================

.. _ATLAS: https://github.com/metagenome-atlas/atlas

.. contents::
    :backlinks: none

This is a alpha version
-----------------------


In this alpha version of atlas the repoprts don't work properly.
I hope they will be available soon.

Most information can be found either in the stats folder or in the input files (see below).


I protected all input files so they don't get deleted and you will be able to re-create the reports onece the beta version is out.

Input files for this report:
{input_files}



"""

    atlas_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    report(
        report_str,
        report_out,
        stylesheet=os.path.join(atlas_dir, "report", "report.css"),
    )


if __name__ == "__main__":

    input_files = []
    for file in snakemake.input:

        logger.info(f"protect file {file}")
        current_file_status = stat.S_IMODE(os.lstat(file).st_mode)
        os.chmod(file, current_file_status & ~stat.S_IEXEC)
        input_files.append(file)

    main(report_out=snakemake.output.report, input_files=input_files)
