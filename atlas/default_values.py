# global defaults
MEM = 80
JAVA_MEM_FRACTION = 0.85
PREALLOCATE_RAM = "t"
PREPROCESS_ADAPTER_MIN_K = 8
PREPROCESS_KMER_TRIM = "r"
PREPROCESS_MINIMUM_BASE_QUALITY = 10
PREPROCESS_ALLOWABLE_KMER_MISMATCHES = 1
PREPROCESS_REFERENCE_KMER_MATCH_LENGTH = 27
QTRIM = "rl"
PREPROCESS_MINIMUM_PASSING_READ_LENGTH = 51
PREPROCESS_MINIMUM_BASE_FREQUENCY = 0.05
PREPROCESS_MAX_NS = -1

MERGING_FLAGS = "ecct iterations=1"
MERGING_EXTEND2 = 50
MERGING_K = 62

CONTAMINANT_MAX_INDEL = 20
CONTAMINANT_MIN_RATIO = 0.65
CONTAMINANT_MINIMUM_HITS = 1
CONTAMINANT_AMBIGUOUS = "best"
CONTAMINANT_KMER_LENGTH = 13

DUPLICATES_ONLY_OPTICAL = False
DUPLICATES_ALLOW_SUBSTITUTIONS = 2

NORMALIZATION_KMER_LENGTH = 21

# almost no filtering unless grossly over-represented
NORMALIZATION_TARGET_DEPTH = 1000  # 500
# allow very low represented kmers to remain
NORMALIZATION_MINIMUM_KMERS = 3  # 15

ASSEMBLY_MEMORY = 250
ASSEMBLY_THREADS = 8
MEGAHIT_MIN_COUNT = 2
MEGAHIT_K_MIN = 21
MEGAHIT_K_MAX = 121
MEGAHIT_K_STEP = 20
MEGAHIT_MERGE_LEVEL = "20,0.98"
MEGAHIT_PRUNE_LEVEL = 2
MEGAHIT_LOW_LOCAL_RATIO = 0.2
SPADES_K = "auto"
# basically no filtering after assembly
PREFILTER_MINIMUM_CONTIG_LENGTH = 200  # 500
# this is bumped up slightly to filter non-merged R1 and R2 sequences
MINIMUM_CONTIG_LENGTH = 300  # 2200

# leave all contigs
MINIMUM_AVERAGE_COVERAGE = 1  # 5
MINIMUM_PERCENT_COVERED_BASES = 20  # 40
MINIMUM_MAPPED_READS = 0
CONTIG_TRIM_BP = 0  # 100

# bases
MINIMUM_REGION_OVERLAP = 1
FEATURE_COUNTS_ALLOW_OVERLAP = True
MAXIMUM_COUNTED_MAP_SITES = 10
# default bbmap
CONTIG_MIN_ID = 0.76
CONTIG_MAP_PAIRED_ONLY = True
CONTIG_MAX_DISTANCE_BETWEEN_PAIRS = 1000
# only best
CONTIG_COUNT_MULTI_MAPPED_READS = False
PROKKA_KINGDOM = "Bacteria"

MAXBIN_MAX_ITERATION = 50
MAXBIN_MIN_CONTIG_LENGTH = 1000
MAXBIN_PROB_THRESHOLD = 0.9

DIAMOND_TOP_SEQS = 2
DIAMOND_E_VALUE = 0.000001
DIAMOND_MIN_IDENTITY = 50
DIAMOND_QUERY_COVERAGE = 60
DIAMOND_GAP_OPEN = 11
DIAMOND_GAP_EXTEND = 1
DIAMOND_BLOCK_SIZE = 2
DIAMOND_INDEX_CHUNKS = 4

SUMMARY_METHOD = "lca"
AGGREGATION_METHOD = "lca-majority"
MAJORITY_THRESHOLD = 0.51
MIN_BITSCORE = 0
MIN_LENGTH = 20
MAX_HITS = 100


import multiprocessing
import os
import sys
import tempfile

EGGNOG_HEADER = [
    "Query",
    "Seed",
    "Seed_evalue",
    "Seed_Score",
    "eggNOG",
    "max_annot_lvl",
    "COG_cat",
    "Description",
    "Name",
    "GO_terms",
    "EC",
    "KO",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
]


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
    config["preprocess_adapter_min_k"] = PREPROCESS_ADAPTER_MIN_K
    config["preprocess_minimum_base_quality"] = PREPROCESS_MINIMUM_BASE_QUALITY
    config[
        "preprocess_allowable_kmer_mismatches"
    ] = PREPROCESS_ALLOWABLE_KMER_MISMATCHES
    config[
        "preprocess_reference_kmer_match_length"
    ] = PREPROCESS_REFERENCE_KMER_MATCH_LENGTH
    config[
        "preprocess_minimum_passing_read_length"
    ] = PREPROCESS_MINIMUM_PASSING_READ_LENGTH
    config["preprocess_minimum_base_frequency"] = PREPROCESS_MINIMUM_BASE_FREQUENCY

    config["deduplicate"] = True
    config["error_correction_overlapping_pairs"] = True

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
    config["error_correction_minprob"] = 0.5  # for memory issues only , I think

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
    config["prefilter_minimum_contig_length"] = PREFILTER_MINIMUM_CONTIG_LENGTH
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

    config["maximum_counted_map_sites"] = MAXIMUM_COUNTED_MAP_SITES

    # gene cluster
    config["genecatalog"] = {
        "source": "genomes",
        "clustermethod": "linclust",
        "minlength": 100,
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

    config["cobining_min_contig_length"] = 2000
    config["cobining_min_bin_size"] = 200 * 1000
    config["semibin_options"] = " --recluster --max-node 1 --max-edges 200 "
    config["cobinning_separator"] = ":"

    config["annotations"] = ["gtdb_taxonomy", "checkm_taxonomy", "gtdb_tree"]
    config["rename_mags_contigs"] = True

    config["genome_dereplication"] = dict(
        filter=dict(noFilter=False, length=5000, completeness=50, contamination=10),
        score=dict(
            completeness=1,
            contamination=5,
            strain_heterogeneity=0,  # not in table
            N50=0.5,
            length=0,
        ),
        ANI=0.95,
        overlap=0.6,
        sketch_size=5000,
        opt_parameters="",
    )

    config["cat_range"] = 5
    config["cat_fraction"] = 0.3

    config["runtime"] = {"default": 5, "assembly": 24, "long": 12, "simple_job": 1}

    # config["diamond_run_mode"] = "fast"
    # config["diamond_top_seqs"] = DIAMOND_TOP_SEQS
    # config["diamond_e_value"] = DIAMOND_E_VALUE
    # config["diamond_min_identity"] = DIAMOND_MIN_IDENTITY
    # config["diamond_query_coverage"] = DIAMOND_QUERY_COVERAGE
    # config["diamond_gap_open"] = DIAMOND_GAP_OPEN
    # config["diamond_gap_extend"] = DIAMOND_GAP_EXTEND
    # config["diamond_block_size"] = DIAMOND_BLOCK_SIZE
    # config["diamond_index_chunks"] = DIAMOND_INDEX_CHUNKS
    # config["summary_method"] = SUMMARY_METHOD
    # config["aggregation_method"] = AGGREGATION_METHOD
    # config["majority_threshold"] = MAJORITY_THRESHOLD

    return config
