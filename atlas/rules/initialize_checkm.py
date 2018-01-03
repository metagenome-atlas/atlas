#!/usr/bin/env python
# coding=utf-8
"""
Initialize CheckM.
"""

import argparse
import logging
import os
import sys
from subprocess import Popen, PIPE


def run_popen(cmd, response, stderr=None):
    # $ checkm data setRoot
    # It seems that the CheckM data folder has not been set yet or has been removed. Running: 'checkm data setRoot'.
    # Where should CheckM store it's data?
    # Please specify a location or type 'abort' to stop trying:
    # /people/brow015/checkm
    #
    # Path [/people/brow015/checkm] has been created and you have permission to write to this folder.
    # (re) creating manifest file (please be patient).
    #
    # You can run 'checkm data update' to ensure you have the latest data files.
    #
    #
    # *******************************************************************************
    #  [CheckM - data] Check for database updates. [setRoot]
    # *******************************************************************************
    #
    # Where should CheckM store it's data?
    # Please specify a location or type 'abort' to stop trying:
    # /people/brow015/checkm
    #
    # Path [/people/brow015/checkm] exists and you have permission to write to this folder.
    # (re) creating manifest file (please be patient).
    #
    # You can run 'checkm data update' to ensure you have the latest data files.
    #
    # Data location successfully changed to: /people/brow015/checkm

    # $ checkm data update
    #
    # *******************************************************************************
    #  [CheckM - data] Check for database updates. [update]
    # *******************************************************************************
    #
    # Connecting to ACE server.
    #
    # ****************************************************************
    # 32 new file(s) to be downloaded from source
    # 0 existing file(s) to be updated
    # 1.39 GB will need to be downloaded
    # Confirm you want to download this data
    # Changes *WILL* be permanent
    # Continue? (y,n) : y
    # ****************************************************************
    newline = os.linesep
    p = Popen(cmd, stdin=PIPE, stdout=PIPE, universal_newlines=True, stderr=stderr)
    stdout = p.communicate(input=newline.join(response) if isinstance(response, list) else response)[0]


def main(db_dir, confirmation, log):
    with open(log, "w") as errlog:
        logging.info("Updating CheckM's data directory.")
        # if checkm has never been run before it prompts for a directory
        # prior to running anything. And then it runs the update code and
        # requires a second input (of the same thing)
        run_popen(["checkm", "data", "setRoot"], [db_dir, db_dir], stderr=errlog)
        # logging.info("Downloading CheckM reference data.")
        # run_popen(["checkm", "data", "update"], ["y", "y"], stderr=errlog)

    # when re-activating a conda env, reset the .dmanifest directory and download
    with open(confirmation, "w") as fh:
        logging.info("CheckM has been successfully initialized.")
        fh.write("CheckM successfully initialized.\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument("dbdir", help="path to database directory")
    p.add_argument("confirmation", help="path to touched file confirming completion")
    p.add_argument("log", help="stderr output of checkm setup")
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, datefmt="%Y-%m-%d %H:%M",
        format="[%(asctime)s] %(message)s")
    main(args.dbdir, args.confirmation, args.log)
