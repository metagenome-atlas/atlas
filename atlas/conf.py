import logging
import multiprocessing
import os
import sys
import tempfile
#import yaml
from ruamel.yaml import YAML
from collections import OrderedDict
from snakemake.io import load_configfile
# default globals
from atlas.default_values import *


writer= YAML()


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



def make_default_config():
    """ generates a dict with all the default values, if they exist
    """

    conf = {}
    conf["tmpdir"] = tempfile.gettempdir()
    conf["threads"] = multiprocessing.cpu_count()
    conf["java_mem"] = JAVA_MEM
    conf["preprocess_adapter_min_k"] = PREPROCESS_ADAPTER_MIN_K
    conf["preprocess_minimum_base_quality"] = PREPROCESS_MINIMUM_BASE_QUALITY
    conf["preprocess_allowable_kmer_mismatches"] = PREPROCESS_ALLOWABLE_KMER_MISMATCHES
    conf["preprocess_reference_kmer_match_length"] = PREPROCESS_REFERENCE_KMER_MATCH_LENGTH
    conf["preprocess_minimum_passing_read_length"] = PREPROCESS_MINIMUM_PASSING_READ_LENGTH
    conf["preprocess_minimum_base_frequency"] = PREPROCESS_MINIMUM_BASE_FREQUENCY

    conf["deduplicate"] = True
    conf["error_correction_overlapping_pairs"] = True
    conf["merge_pairs_before_assembly"] = True

    conf["contaminant_max_indel"] = CONTAMINANT_MAX_INDEL
    conf["contaminant_min_ratio"] = CONTAMINANT_MIN_RATIO
    conf["contaminant_kmer_length"] = CONTAMINANT_KMER_LENGTH
    conf["contaminant_minimum_hits"] = CONTAMINANT_MINIMUM_HITS
    conf["contaminant_ambiguous"] = CONTAMINANT_AMBIGUOUS

    conf["duplicates_only_optical"] = DUPLICATES_ONLY_OPTICAL
    conf["duplicates_allow_substitutions"] = DUPLICATES_ALLOW_SUBSTITUTIONS

    conf["normalization_kmer_length"] = NORMALIZATION_KMER_LENGTH
    conf["normalization_target_depth"] = NORMALIZATION_TARGET_DEPTH
    conf["normalization_minimum_kmers"] = NORMALIZATION_MINIMUM_KMERS

    conf["merging_k"] = MERGING_K
    conf["merging_extend2"] = MERGING_EXTEND2
    conf["merging_flags"]  = MERGING_FLAGS

    conf["assembler"] = 'megahit'
    conf["assembly_memory"] = ASSEMBLY_MEMORY
    conf["assembly_threads"] = ASSEMBLY_THREADS
    conf["megahit_min_count"] = MEGAHIT_MIN_COUNT
    conf["megahit_k_min"] = MEGAHIT_K_MIN
    conf["megahit_k_max"] = MEGAHIT_K_MAX
    conf["megahit_k_step"] = MEGAHIT_K_STEP
    conf["megahit_merge_level"] = MEGAHIT_MERGE_LEVEL
    conf["megahit_prune_level"] = MEGAHIT_PRUNE_LEVEL
    conf["megahit_low_local_ratio"] = MEGAHIT_LOW_LOCAL_RATIO
    conf["minimum_contig_length"] = MINIMUM_CONTIG_LENGTH
    conf["prefilter_minimum_contig_length"] = PREFILTER_MINIMUM_CONTIG_LENGTH
    conf["spades_k"] = SPADES_K
    conf["minimum_average_coverage"] = MINIMUM_AVERAGE_COVERAGE
    conf["minimum_percent_covered_bases"] = MINIMUM_PERCENT_COVERED_BASES
    conf["minimum_mapped_reads"] = MINIMUM_MAPPED_READS
    conf["contig_trim_bp"] = CONTIG_TRIM_BP

    conf["translation_table"] = 11

    # map reads to contigs and to genes
    conf["minimum_region_overlap"] = MINIMUM_REGION_OVERLAP
    conf["feature_counts_allow_overlap"] = FEATURE_COUNTS_ALLOW_OVERLAP
    conf["contig_count_multi_mapped_reads"] = CONTIG_COUNT_MULTI_MAPPED_READS
    conf["contig_min_id"] = CONTIG_MIN_ID
    conf["contig_map_paired_only"] = CONTIG_MAP_PAIRED_ONLY
    conf["contig_max_distance_between_pairs"] = CONTIG_MAX_DISTANCE_BETWEEN_PAIRS


    conf["maximum_counted_map_sites"] = MAXIMUM_COUNTED_MAP_SITES
    conf["perform_genome_binning"] = True
    conf["maxbin_max_iteration"] = MAXBIN_MAX_ITERATION
    conf["maxbin_min_contig_length"] = MAXBIN_MIN_CONTIG_LENGTH
    conf["maxbin_prob_threshold"] = MAXBIN_PROB_THRESHOLD

    conf["diamond_run_mode"] = "fast"
    conf["diamond_top_seqs"] = DIAMOND_TOP_SEQS
    conf["diamond_e_value"] = DIAMOND_E_VALUE
    conf["diamond_min_identity"] = DIAMOND_MIN_IDENTITY
    conf["diamond_query_coverage"] = DIAMOND_QUERY_COVERAGE
    conf["diamond_gap_open"] = DIAMOND_GAP_OPEN
    conf["diamond_gap_extend"] = DIAMOND_GAP_EXTEND
    conf["diamond_block_size"] = DIAMOND_BLOCK_SIZE
    conf["diamond_index_chunks"] = DIAMOND_INDEX_CHUNKS
    conf["summary_method"] = SUMMARY_METHOD
    conf["aggregation_method"] = AGGREGATION_METHOD
    conf["majority_threshold"] = MAJORITY_THRESHOLD

    return conf

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

    template_conf_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),"template_config.yaml")

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
