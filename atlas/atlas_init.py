import os, sys
from .color_logger import logger
import multiprocessing
import tempfile
from snakemake import utils
from snakemake.io import load_configfile
import pandas as pd
import numpy as np
from collections import defaultdict
import click

from .make_config import make_config, validate_config
from .create_sample_table import get_samples_from_fastq,simplify_sample_names
from .sample_table import validate_sample_table, load_sample_table,ADDITIONAL_SAMPLEFILE_HEADERS

# default globals
ADAPTERS = "adapters.fa"
RRNA = "silva_rfam_all_rRNAs.fa"
PHIX = "phiX174_virus.fa"
SRA_READ_PATH = "SRAreads"





def prepare_sample_table_for_atlas(
    sample_table, reads_are_QC=False, outfile="samples.tsv"
):

    """
    Write the file `samples.tsv` and complete the sample names and paths for all
    files in `path`.
    Args:
            path_to_fastq (str): fastq/fasta data directory
    """



    if os.path.exists(outfile):
        logger.error(
            f"Output file {outfile} already exists I don't dare to overwrite it."
        )
        exit(1)

    simplify_sample_names(sample_table)

    columns = sample_table.columns  # R1 and R2 or only R1 , who knows

    # Test if paired end
    if "R2" in columns:
        fractions=["R1","R2"]
    else:
        sample_table.rename(columns = {"R1":"se"}, inplace=True)
        fractions=["se"]

    # Add prefix to fractions depending if qc or not
    if reads_are_QC:
        prefix= "Reads_QC_"
    else:
        prefix= "Reads_raw_"

    sample_table.rename(columns = {f:f"{prefix}{f}" for f in fractions},inplace=True)

    # Add BinGroup and additional empty headers
    Headers = ADDITIONAL_SAMPLEFILE_HEADERS
    for h in Headers:
        sample_table[h] = np.nan

    sample_table["BinGroup"] = sample_table.index

    validate_sample_table(sample_table)




    sample_table.to_csv(outfile, sep="\t")






@click.command(
    "init",
    short_help="prepare configuration file and sample table for atlas run",
)
@click.argument("path_to_fastq", type=click.Path(readable=True))
@click.option(
    "-d",
    "--db-dir",
    default=os.path.join(os.path.realpath("."), "databases"),
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    show_default=True,
    help="location to store databases (need ~50GB)",
)
@click.option(
    "-w",
    "--working-dir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="location to run atlas",
    default=".",
)
@click.option(
    "--assembler",
    default="spades",
    type=click.Choice(["megahit", "spades"]),
    show_default=True,
    help="assembler",
)
@click.option(
    "--data-type",
    default="metagenome",
    type=click.Choice(["metagenome", "metatranscriptome"]),
    show_default=True,
    help="sample data type",
)
@click.option(
    "--interleaved-fastq",
    is_flag=True,
    default=False,
    help="fastq files are paired-end in one files (interleaved)",
)
@click.option(
    "--threads",
    default=8,
    type=int,
    help="number of threads to use per multi-threaded job",
)
@click.option(
    "--skip-qc",
    is_flag=True,
    help="Skip QC, if reads are already pre-processed",
)
def run_init(
    path_to_fastq,
    db_dir,
    working_dir,
    assembler,
    data_type,
    interleaved_fastq,
    threads,
    skip_qc=False,
):
    """Write the file CONFIG and complete the sample names and paths for all
    FASTQ files in PATH.

    PATH is traversed recursively and adds any file with '.fastq' or '.fq' in
    the file name with the file name minus extension as the sample ID.
    """

    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    config = os.path.join(working_dir, "config.yaml")
    if not os.path.exists(db_dir):
        os.makedirs(db_dir)
    sample_file = os.path.join(working_dir, "samples.tsv")

    make_config(db_dir, threads, assembler, data_type, interleaved_fastq, config)
    sample_table = get_samples_from_fastq(path_to_fastq)
    prepare_sample_table_for_atlas(
        sample_table, reads_are_QC=skip_qc, outfile=sample_file
    )


@click.command(
    "init-public",
    short_help="Prepare atlas run from public data from SRA",
    help="prepare configuration file and sample table for atlas run"
    "based on public data from SRA\n"
    "Supply a set of SRA run ids to the command:"
    "SRR4305427 ERR1190946\n\n"
    "Reads are automatically downloaded but not stored on your machine. Beware that sometimes multiple runs go into one sample.",
)
@click.argument("identifiers", nargs=-1)
@click.option(
    "-d",
    "--db-dir",
    default=os.path.join(os.path.realpath("."), "databases"),
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    show_default=True,
    help="location to store databases (need ~50GB)",
)
@click.option(
    "-w",
    "--working-dir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="location to run atlas",
    default=".",
)
@click.option(
    "--skip-qc",
    is_flag=True,
    help="Skip QC, if reads are already pre-processed",
)
# @click.option(
#     "--single-end",
#     is_flag=True,
#     help="Your reads are single end",
# )
def run_init_sra(identifiers, db_dir, working_dir, skip_qc=False, single_end=False):
    """Write the file CONFIG and complete the sample.tsv"""

    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    config = os.path.join(working_dir, "config.yaml")
    if not os.path.exists(db_dir):
        os.makedirs(db_dir)
    sample_file = os.path.join(working_dir, "samples.tsv")

    make_config(db_dir, config=config)
    sample_table = pd.DataFrame(index=identifiers)

    if single_end:
        sample_table["R1"] = sample_table.index.map(
            lambda s: os.path.join(SRA_READ_PATH, f"{s}.fastq.gz")
        )
    else:

        sample_table["R1"] = sample_table.index.map(
            lambda s: os.path.join(SRA_READ_PATH, f"{s}_1.fastq.gz")
        )
        sample_table["R2"] = sample_table.index.map(
            lambda s: os.path.join(SRA_READ_PATH, f"{s}_2.fastq.gz")
        )
        prepare_sample_table_for_atlas(
            sample_table, reads_are_QC=skip_qc, outfile=sample_file
        )
