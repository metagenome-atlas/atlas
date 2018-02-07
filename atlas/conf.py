import logging
import multiprocessing
import os
import sys
import tempfile
from ruamel.yaml import YAML
from snakemake.io import load_configfile
# default globals


def get_sample_files(path, data_type):
    samples = dict()
    seen = set()
    for dir_name, sub_dirs, files in os.walk(path):
        for fname in files:

            if ".fastq" in fname or ".fq" in fname:

                sample_id = fname.partition(".fastq")[0]
                if ".fq" in sample_id:
                    sample_id = fname.partition(".fq")[0]

                sample_id = sample_id.replace("_R1", "").replace("_r1", "").replace("_R2", "").replace("_r2", "")
                sample_id = sample_id.replace("_", "-").replace(" ", "-")

                fq_path = os.path.join(dir_name, fname)
                fastq_paths = [fq_path]

                if fq_path in seen: continue

                if "_R1" in fname or "_r1" in fname:
                    r2_path = os.path.join(dir_name, fname.replace("_R1", "_R2").replace("_r1", "_r2"))
                    if not r2_path == fq_path:
                        seen.add(r2_path)
                        fastq_paths.append(r2_path)

                if "_R2" in fname or "_r2" in fname:
                    r1_path = os.path.join(dir_name, fname.replace("_R2", "_R1").replace("_r2", "_r1"))
                    if not r1_path == fq_path:
                        seen.add(r1_path)
                        fastq_paths.insert(0, r1_path)

                if sample_id in samples:
                    logging.warn("Duplicate sample %s was found after renaming; skipping..." % sample_id)
                    continue
                samples[sample_id] = {'fastq': fastq_paths, 'type': data_type}
    return samples


def make_config(config, path, data_type, database_dir, threads, assembler):
    """
    Reads template config file with comments from ./template_config.yaml
    updates it by the parameters provided.
    Write the file `config` and complete the sample names and paths for all
    files in `path`.

    Args:
        config (str): output file path for yaml
        path (str): fastq/fasta data directory
        data_type (str): this is either metagenome or metatranscriptome
        database_dir (str): location of downloaded databases
        threads (int): number of threads per node to utilize
        assembler (str): either spades or megahit
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

    conf["samples"] = samples
    conf["tmpdir"] = tempfile.gettempdir()
    conf["threads"] = multiprocessing.cpu_count() if not threads else threads
    conf["preprocess_adapters"] = os.path.join(database_dir, ADAPTERS)
    conf["contaminant_references"] = {"rRNA":os.path.join(database_dir, RRNA),
                                      "PhiX":os.path.join(database_dir, PHIX)}

    conf["assembler"] = assembler

    conf["refseq_namemap"] = os.path.join(database_dir, "refseq.db")
    conf["refseq_tree"] = os.path.join(database_dir, "refseq.tree")
    conf["diamond_db"] = os.path.join(database_dir, "refseq.dmnd")

    with open(config, "w") as f:
        yaml.dump(conf, f)
    logging.info("Configuration file written to %s" % config)


def log_exception(msg):
    logging.critical(msg)
    logging.info("Documentation is available at: https://pnnl-atlas.readthedocs.io")
    logging.info("Issues can be raised at: https://github.com/pnnl/atlas/issues")
    sys.exit(1)


def validate_sample_defs(config):
    if "samples" not in config.keys():
        log_exception("'samples' are not defined in the configuration")
    if len(config["samples"]) == 0:
        log_exception("No samples are defined under 'samples'")
    for sample in config["samples"]:
        if "fastq" not in config["samples"][sample]:
            log_exception("'fastq' must be defined per sample")
        # single- or paired-end and user added appropriately as list
        if isinstance(config["samples"][sample]["fastq"], list):
            for fq in config["samples"][sample]["fastq"]:
                if not os.path.exists(fq):
                    log_exception("%s does not exist" % fq)
        # single-end and user defined as string
        else:
            if not os.path.exists(config["samples"][sample]["fastq"]):
                log_exception("%s does not exist" % config["samples"][sample]["fastq"])


def validate_config(config):
    conf = load_configfile(config)
    validate_sample_defs(conf)
    # could later add more validation steps
