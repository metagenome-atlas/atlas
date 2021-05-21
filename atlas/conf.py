import logging
import multiprocessing
import os
import sys
import tempfile
from snakemake import utils
from snakemake.io import load_configfile
import pandas as pd
import numpy as np
from collections import defaultdict
import click

sys.path.append(os.path.dirname(__file__))
from default_values import make_default_config

# default globals
ADAPTERS = "adapters.fa"
RRNA = "silva_rfam_all_rRNAs.fa"
PHIX = "phiX174_virus.fa"
ADDITIONAL_SAMPLEFILE_HEADERS = []  # ,'Contigs']


def get_samples_from_fastq(path):
    """
    creates table sampleID R1 R2 with the absolute paths of fastq files in a given folder
    """
    samples = defaultdict(dict)
    seen = set()
    for dir_name, sub_dirs, files in os.walk(os.path.abspath(path)):
        for fname in files:

            if ".fastq" in fname or ".fq" in fname:

                sample_id = fname.split(".fastq")[0].split(".fq")[0]

                sample_id = (
                    sample_id.replace("_R1", "")
                    .replace("_r1", "")
                    .replace("_R2", "")
                    .replace("_r2", "")
                )
                sample_id = sample_id.replace("_", "-").replace(" ", "-")

                fq_path = os.path.join(dir_name, fname)

                if fq_path in seen:
                    continue

                if "_R2" in fname or "_r2" in fname:

                    if "R2" in samples[sample_id]:
                        logging.error(
                            f"Duplicate sample {sample_id} was found after renaming; skipping... \n Samples: \n{samples}"
                        )

                    samples[sample_id]["R2"] = fq_path
                else:
                    if "R1" in samples[sample_id]:
                        logging.error(
                            f"Duplicate sample {sample_id} was found after renaming; skipping... \n Samples: \n{samples}"
                        )

                    samples[sample_id]["R1"] = fq_path

    samples = pd.DataFrame(samples).T

    if samples.isnull().any().any():
        logging.error(f"Missing files:\n\n {samples}")
        exit(1)

    if samples.shape[0] == 0:
        logging.error(
            f"No files found in {path}\n"
            "I'm looking for files with .fq or .fastq extension. "
        )
        exit(1)

    return samples


def validate_sample_table(sampleTable):

    Expected_Headers = ["BinGroup"] + ADDITIONAL_SAMPLEFILE_HEADERS
    for h in Expected_Headers:
        if not (h in sampleTable.columns):
            logging.error(f"expect '{h}' to be found in samples.tsv")
            exit(1)
        elif sampleTable[h].isnull().any():
            logging.error(f"Found empty values in the sample table column '{h}'")
            exit(1)

    if not sampleTable.index.is_unique:
        duplicated_samples = ", ".join(D.index.duplicated())
        logging.error(
            f"Expect Samples to be unique. Found {duplicated_samples} more than once"
        )
        exit(1)


def prepare_sample_table(path_to_fastq, reads_are_QC=False, outfile="samples.tsv"):
    """
    Write the file `samples.tsv` and complete the sample names and paths for all
    files in `path`.
    Args:
            path_to_fastq (str): fastq/fasta data directory
    """

    samples = get_samples_from_fastq(path_to_fastq)

    columns = samples.columns  # R1 and R2 or only R1 , who knows

    if "R2" not in columns:
        assert len(columns) == 1, "expect columns to be only ['R1']"
        columns = ["se"]

    if reads_are_QC:
        samples.columns = ["Reads_QC_" + c for c in columns]
    else:
        samples.columns = ["Reads_raw_" + c for c in columns]

    Headers = ADDITIONAL_SAMPLEFILE_HEADERS

    for h in Headers:
        samples[h] = np.nan

    samples["BinGroup"] = samples.index

    validate_sample_table(samples)

    logging.info("Found %d samples under %s" % (len(samples), path_to_fastq))
    if os.path.exists(outfile):
        logging.error(
            f"Output file {outfile} already exists I don't dare to overwrite it."
        )
        exit(1)
    else:
        samples.to_csv(outfile, sep="\t")


def load_sample_table(sample_table="samples.tsv"):

    sampleTable = pd.read_csv(sample_table, index_col=0, sep="\t")
    validate_sample_table(sampleTable)
    return sampleTable


def make_config(
    database_dir,
    threads,
    assembler,
    data_type="metagenome",
    interleaved_fastq=False,
    config="config.yaml",
):
    """
    Reads template config file with comments from ./template_config.yaml
    updates it by the parameters provided.

    Args:
        config (str): output file path for yaml
        database_dir (str): location of downloaded databases
        threads (int): number of threads per node to utilize
        assembler (str): either spades or megahit
        data_type (str): this is either metagenome or metatranscriptome
    """

    from ruamel.yaml import YAML  # used for yaml reading with comments

    yaml = YAML()
    yaml.version = (1, 1)
    yaml.default_flow_style = False

    template_conf_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "template_config.yaml"
    )

    with open(template_conf_file) as template_config:
        conf = yaml.load(template_config)

    conf["tmpdir"] = tempfile.gettempdir()
    conf["threads"] = multiprocessing.cpu_count() if not threads else threads
    conf["preprocess_adapters"] = os.path.join(database_dir, "adapters.fa")
    conf["contaminant_references"] = {
        "PhiX": os.path.join(database_dir, "phiX174_virus.fa")
    }

    if data_type == "metatranscriptome":
        conf["contaminant_references"]["rRNA"] = (
            os.path.join(database_dir, "silva_rfam_all_rRNAs.fa"),
        )

    conf["data_type"] = data_type
    conf["interleaved_fastqs"] = interleaved_fastq

    conf["assembler"] = assembler
    conf["database_dir"] = database_dir
    # conf["refseq_namemap"] = os.path.join(database_dir, "refseq.db")
    # conf["refseq_tree"] = os.path.join(database_dir, "refseq.tree")
    # conf["diamond_db"] = os.path.join(database_dir, "refseq.dmnd")

    if os.path.exists(config):
        logging.warning(
            f"Config file {config} already exists, I didn't dare to overwrite it. continue..."
        )
    else:

        with open(config, "w") as f:
            yaml.dump(conf, f)
        logging.info(
            "Configuration file written to %s\n"
            "You may want to edit it using any text editor." % config
        )


def validate_config(config, workflow):
    conf = load_configfile(config)


#    validate_sample_defs(conf, workflow)
# could later add more validation steps


def update_config(config):
    """
    Populates config file with default config values.
    And made changes if necessary.

    """

    # in old version java_mem was used, new is mem
    if ("java_mem" in config) and (not ("mem" in config)):
        config["mem"] = config["java_mem"]

    # get default values and update them with values specified in config file
    default_config = make_default_config()
    utils.update_config(default_config, config)

    return default_config


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
    prepare_sample_table(path_to_fastq, reads_are_QC=skip_qc, outfile=sample_file)
