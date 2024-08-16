#!/usr/bin/env python
import os
import requests
import csv, io
import traceback
import xmltodict

import pandas as pd

import logging


logger = logging.Logger("public-init")


from .kingfisher import SraMetadata

SRA = SraMetadata()


def get_table_from_accessions(accessions, output_file="SRA_run_table.csv") -> None:

    accessions = set(accessions)

    bioproject_accessions = [
        identifier
        for identifier in accessions
        if identifier.startswith("PRJ") or identifier[1:3] == "RP"
    ]

    run_accessions = [
        identifier for identifier in accessions if identifier[1:3] == "RR"
    ]

    not_searchable = list(accessions - set(run_accessions) - set(bioproject_accessions))
    n_searchable = len(bioproject_accessions) + len(run_accessions)

    if n_searchable < len(accessions):

        error_text = "I can only search for run-accessions, e.g. ERR1739691 "
        " or projects e.g. PRJNA621514 or SRP260223\n"

        if n_searchable > 0:
            error_text += f"I can search for {n_searchable} accessions,\n but "

        error_text += (
            f"I cannot search for {len(not_searchable)} accessions: "
            + ",".join(not_searchable[0 : min(len(not_searchable), 4)])
        )

        raise Exception(error_text)

    # get sra for project

    if len(bioproject_accessions) > 0:

        run_accessions_for_bioprojects = SRA.fetch_runs_from_bioprojects(
            bioproject_accessions
        )

        logger.info(
            f"Found {len(run_accessions_for_bioprojects)} runs for {len(bioproject_accessions)} bioprojects"
        )

        overlap = set(run_accessions) & set(run_accessions_for_bioprojects)

        if len(overlap) > 0:
            logger.warning(
                f"Found {len(overlap)} runs that are in both the bioprojects and the run accessions, downloading only once"
            )
            run_accessions_for_bioprojects = list(
                set(run_accessions_for_bioprojects) - overlap
            )

        run_accessions = run_accessions + run_accessions_for_bioprojects

    # Fetch runs
    logger.info(f"Fetch info for {len(run_accessions)} runs")

    tab = SRA.efetch_sra_from_accessions(run_accessions)

    tab.set_index("run", inplace=True)

    tab.to_csv(output_file, index=True)

    logger.info(f"Table saved to {output_file}")


def parse_arguments_from_terminal():
    ## Comand line interface
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "accessions",
        type=str,
        help="A list of space separated SRA identifiers runs or bioprojects e.g. PRJNA63661 SRP046206 SRR1759594",
        nargs="+",
    )
    parser.add_argument(
        "-o",
        "--outprefix",
        type=str,
        help="Output prefix for the sra table",
        default="SRA",
    )
    parser.add_argument(
        "-f",
        "--overwrite",
        action="store_true",
        help="Overwrite table if it already exists",
    )

    # parse arguments from terminal
    args = parser.parse_args()

    # if output prefix is only a folder add 'SRA
    outprefix = args.outprefix
    if outprefix.endswith(os.path.sep):
        outprefix += "SRA"

    # Define output_file
    output_file = outprefix + "_runtable.csv"

    if os.path.exists(output_file) and not args.overwrite:
        logger.error(
            f"File `{output_file}` already exists, use --overwrite to over write it or choose another prefix"
        )
        exit(1)

    accessions = args.accessions
    # If one string transform
    if type(accessions) == str:
        accessions = [accessions]

    return accessions, output_file


if __name__ == "__main__":
    get_table_from_accessions(*parse_arguments_from_terminal())
