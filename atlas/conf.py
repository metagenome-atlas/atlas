import logging
import multiprocessing
import os
import sys
import tempfile
from ruamel.yaml import YAML
from snakemake.io import load_configfile
import logging
import pandas as pd
import numpy as np
from collections import defaultdict
# default globals
ADAPTERS = "adapters.fa"
RRNA = "silva_rfam_all_rRNAs.fa"
PHIX = "phiX174_virus.fa"
ADDITIONAL_SAMPLEFILE_HEADERS=['Source_group','Contigs']

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

                sample_id = sample_id.replace("_R1", "").replace("_r1", "").replace("_R2", "").replace("_r2", "")
                sample_id = sample_id.replace("_", "-").replace(" ", "-")

                fq_path = os.path.join(dir_name, fname)

                if fq_path in seen: continue

                if "_R2" in fname or "_r2" in fname:

                    if 'R2' in samples[sample_id]:
                        logging.error(f"Duplicate sample {sample_id} was found after renaming; skipping... \n Samples: \n{samples}")

                    samples[sample_id]['R2'] = fq_path
                else:
                    if 'R1' in samples[sample_id]:
                        logging.error(f"Duplicate sample {sample_id} was found after renaming; skipping... \n Samples: \n{samples}")

                    samples[sample_id]['R1'] = fq_path


    samples= pd.DataFrame(samples).T

    if samples.isna().any().any():
        logging.error(f"Missing files:\n\n {samples}")

    return samples

def prepare_sample_table(path_to_fastq,reads_are_QC=False,outfile='samples.tsv'):
    """
    Write the file `samples.tsv` and complete the sample names and paths for all
    files in `path`.
    Args:
            path_to_fastq (str): fastq/fasta data directory
    """

    samples = get_samples_from_fastq(path_to_fastq)

    if reads_are_QC:
        samples.columns= 'Reads_QC_'+samples.columns
    else:
        samples.columns= 'Reads_raw_'+samples.columns
        ADDITIONAL_SAMPLEFILE_HEADERS = list('Reads_QC_'+samples.columns) + ADDITIONAL_SAMPLEFILE_HEADERS

    for h in ADDITIONAL_HEADERS:
        samples[h]=np.nan

    if os.path.exists(outfile):
        logging.error(f"Output file {outfile} already exists I don't dare to overwrite it.")
    else:
        samples.to_csv(outfile,sep='\t')



def make_config(database_dir, threads, assembler, data_type='metagenome',config='config.yaml'):
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
    config = os.path.realpath(os.path.expanduser(config))
    os.makedirs(os.path.dirname(config), exist_ok=True)

    path = os.path.realpath(os.path.expanduser(path))
    database_dir = os.path.realpath(os.path.expanduser(database_dir))

    yaml = YAML()
    yaml.version = (1, 1)
    yaml.default_flow_style = False

    template_conf_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      "template_config.yaml")

    with open(template_conf_file) as template_config:
        conf = yaml.load(template_config)

    samples = get_sample_files(path, data_type)
    logging.info("Found %d samples under %s" % (len(samples), path))

    conf["tmpdir"] = tempfile.gettempdir()
    conf["threads"] = multiprocessing.cpu_count() if not threads else threads
    conf["preprocess_adapters"] = os.path.join(database_dir, "adapters.fa")
    conf["contaminant_references"] = {"rRNA":os.path.join(database_dir, "silva_rfam_all_rRNAs.fa"),
                                      "PhiX":os.path.join(database_dir, "phiX174_virus.fa")}

    #Samples
    conf["data_type"]= data_type


    conf["assembler"] = assembler
    conf["database_dir"] = database_dir
    conf["refseq_namemap"] = os.path.join(database_dir, "refseq.db")
    conf["refseq_tree"] = os.path.join(database_dir, "refseq.tree")
    conf["diamond_db"] = os.path.join(database_dir, "refseq.dmnd")

    with open(config, "w") as f:
        yaml.dump(conf, f)
    logging.info("Configuration file written to %s" % config)


def log_exception(msg):
    logging.critical(msg)
    logging.info("Documentation is available at: https://metagenome-atlas.readthedocs.io")
    logging.info("Issues can be raised at: https://github.com/metagenome-atlas/atlas/issues")
    sys.exit(1)

def validate_config(config, workflow):
    conf = load_configfile(config)
#    validate_sample_defs(conf, workflow)
    # could later add more validation steps
