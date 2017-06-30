import argparse
import os
import shutil
from subprocess import Popen, PIPE


def run_popen(cmd, response):
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
    p = Popen(cmd, stdin=PIPE, stdout=PIPE, universal_newlines=True)
    stdout = p.communicate(input=newline.join(response) if isinstance(response, list) else response)[0]


run_popen(["checkm", "data", "setRoot"], [snakemake.params.database_dir, snakemake.params.database_dir])
run_popen(["checkm", "data", "update"], ["y", "y"])

# when re-activating a conda env, reset the .dmanifest directory and download
with open(snakemake.output.touched_output, "w") as fh:
    fh.write("CheckM successfully initialized.\n")
