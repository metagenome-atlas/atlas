from .default_values import *
from snakemake.utils import update_config as snakemake_update_config
from snakemake.common.configfile import load_configfile
import tempfile
import sys
import os
import multiprocessing
import logging

logger = logging.getLogger(__file__)


def make_default_config():
    """generates a dict with all the default values, if they exist"""

    config = {}
    config["data_type"] = "metagenome"
    config["interleaved_fastqs"] = False
    config["tmpdir"] = tempfile.gettempdir()
    config["threads"] = 6

    config["simplejob_mem"] = 10
    config["simplejob_threads"] = 4
    config[
        "importqc_params"
    ] = "iupacToN=t touppercase=t qout=33 addslash=t trimreaddescription=t"

    config["mem"] = MEM
    config["large_mem"] = 250
    config["large_threads"] = 16
    config["deduplicate"] = True

    config["contaminant_max_indel"] = CONTAMINANT_MAX_INDEL
    config["contaminant_min_ratio"] = CONTAMINANT_MIN_RATIO
    config["contaminant_kmer_length"] = CONTAMINANT_KMER_LENGTH
    config["contaminant_minimum_hits"] = CONTAMINANT_MINIMUM_HITS
    config["contaminant_ambiguous"] = CONTAMINANT_AMBIGUOUS

    config["duplicates_only_optical"] = DUPLICATES_ONLY_OPTICAL
    config["duplicates_allow_substitutions"] = DUPLICATES_ALLOW_SUBSTITUTIONS

    config["normalize_reads_before_assembly"] = False
    config["normalization_kmer_length"] = NORMALIZATION_KMER_LENGTH
    config["normalization_target_depth"] = NORMALIZATION_TARGET_DEPTH
    config["normalization_minimum_kmer_depth"] = 5

    config["error_correction_before_assembly"] = True
    config["error_correction_remove_lowdepth"] = False
    config["error_correction_kmer"] = 31
    config["error_correction_lowdepth_fraction"] = 0.5
    config["error_correction_minimum_kmer_depth"] = 1
    config["error_correction_aggressive"] = False
    # for memory issues only , I think
    config["error_correction_minprob"] = 0.5

    config["merge_pairs_before_assembly"] = True
    config["merging_k"] = MERGING_K
    config["merging_extend2"] = MERGING_EXTEND2
    config["merging_flags"] = MERGING_FLAGS

    config["assembler"] = "spades"
    config["assembly_memory"] = ASSEMBLY_MEMORY
    config["assembly_threads"] = ASSEMBLY_THREADS
    config["megahit_min_count"] = MEGAHIT_MIN_COUNT
    config["megahit_k_min"] = MEGAHIT_K_MIN
    config["megahit_k_max"] = MEGAHIT_K_MAX
    config["megahit_k_step"] = MEGAHIT_K_STEP
    config["megahit_merge_level"] = MEGAHIT_MERGE_LEVEL
    config["megahit_prune_level"] = MEGAHIT_PRUNE_LEVEL
    config["megahit_low_local_ratio"] = MEGAHIT_LOW_LOCAL_RATIO
    config["megahit_preset"] = "default"
    config["minimum_contig_length"] = MINIMUM_CONTIG_LENGTH
    config["spades_k"] = SPADES_K
    config["spades_use_scaffolds"] = True
    config["spades_preset"] = "meta"
    config["spades_extra"] = ""
    config["spades_skip_BayesHammer"] = False
    config["longread_type"] = None
    config["filter_contigs"] = True
    config["minimum_average_coverage"] = MINIMUM_AVERAGE_COVERAGE
    config["minimum_percent_covered_bases"] = MINIMUM_PERCENT_COVERED_BASES
    config["minimum_mapped_reads"] = MINIMUM_MAPPED_READS
    config["contig_trim_bp"] = CONTIG_TRIM_BP

    # config["translation_table"] = 11

    # map reads to contigs and to genes
    config["minimum_region_overlap"] = MINIMUM_REGION_OVERLAP
    config["feature_counts_allow_overlap"] = FEATURE_COUNTS_ALLOW_OVERLAP
    config["contig_count_multi_mapped_reads"] = CONTIG_COUNT_MULTI_MAPPED_READS
    config["contig_min_id"] = CONTIG_MIN_ID
    config["contig_map_paired_only"] = CONTIG_MAP_PAIRED_ONLY
    config["contig_max_distance_between_pairs"] = CONTIG_MAX_DISTANCE_BETWEEN_PAIRS
    config["minimum_map_quality"] = MINIMUM_MAP_QUALITY
    config["maximum_counted_map_sites"] = MAXIMUM_COUNTED_MAP_SITES

    config["bin_quality_asesser"] = "checkm"
    # gene cluster
    config["genecatalog"] = {
        "source": "genomes",
        "clustermethod": "linclust",
        "minlength_nt": 270,
        "minid": 0.9,
        "coverage": 0.9,
        "extra": "",
        "SubsetSize": 500000,
    }

    config["eggNOG_use_virtual_disk"] = False
    config["virtual_disk"] = "/dev/shm"

    # binning
    config["perform_genome_binning"] = True

    config["final_binner"] = "DASTool"
    config["binner"] = ["metabat", "maxbin"]

    config["metabat"] = {"sensitivity": "sensitive", "min_contig_length": 1500}

    config["concoct"] = {
        "Nexpected_clusters": 200,  # important parameter
        "read_length": 100,  # change this parameter !
        "Niterations": 500,
        "min_contig_length": 1000,
    }

    config["maxbin"] = {
        "max_iteration": MAXBIN_MAX_ITERATION,
        "prob_threshold": MAXBIN_PROB_THRESHOLD,
        "min_contig_length": MAXBIN_MIN_CONTIG_LENGTH,
    }

    config["DASTool"] = {
        "search_engine": "diamond",
        "score_threshold": 0.5,
        "duplicate_penalty": 0.6,
        "megabin_penalty": 0.5,
    }

    config["gunc_database"] = "gtdb"

    config["cobining_min_contig_length"] = 2000
    config["cobining_min_bin_size"] = 200 * 1000
    config["cobinning_separator"] = ":"

    config["annotations"] = ["gtdb_taxonomy", "checkm_taxonomy", "gtdb_tree"]
    config["rename_mags_contigs"] = True

    config["runtime"] = {"default": 5, "assembly": 24, "long": 12, "simplejob": 1}

    return config


def make_config(
    database_dir,
    threads=8,
    assembler="spades",
    data_type="metagenome",
    interleaved_fastq=False,
    config="config.yaml",
    binner="vamb",
):
    """
    Reads template config file with comments from ../workflow/config/template_config.yaml
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
    # yaml.version = (1, 1)
    # yaml.default_flow_style = False

    template_conf_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "workflow/../config/template_config.yaml",
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

    conf["final_binner"] = binner

    if os.path.exists(config):
        logger.warning(
            f"Config file {config} already exists, I didn't dare to overwrite it. continue..."
        )
    else:
        with open(config, "w") as f:
            yaml.dump(conf, f)
        logger.info(
            "Configuration file written to %s\n"
            "        You may want to edit it using any text editor." % config
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

    # get default values and update them with values specified in config file
    default_config = make_default_config()
    snakemake_update_config(default_config, config)

    return default_config
