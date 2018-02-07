
# global defaults
JAVA_MEM = 32
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

MERGING_FLAGS = "ecct vstrict"
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
NORMALIZATION_TARGET_DEPTH = 100
NORMALIZATION_MINIMUM_KMERS = 15

ASSEMBLY_MEMORY = 50
ASSEMBLY_THREADS = 8
MEGAHIT_MIN_COUNT = 2
MEGAHIT_K_MIN = 21
MEGAHIT_K_MAX = 121
MEGAHIT_K_STEP = 20
MEGAHIT_MERGE_LEVEL = "20,0.98"
MEGAHIT_PRUNE_LEVEL = 2
MEGAHIT_LOW_LOCAL_RATIO = 0.2
SPADES_K = "auto"
PREFILTER_MINIMUM_CONTIG_LENGTH = 500
MINIMUM_CONTIG_LENGTH = 2200

MINIMUM_AVERAGE_COVERAGE = 5
MINIMUM_PERCENT_COVERED_BASES = 40
MINIMUM_MAPPED_READS = 0
CONTIG_TRIM_BP = 100

# bases
MINIMUM_REGION_OVERLAP = 1
FEATURE_COUNTS_ALLOW_OVERLAP = True
MAXIMUM_COUNTED_MAP_SITES = 10
#default bbmap
CONTIG_MIN_ID = 0.76
CONTIG_MAP_PAIRED_ONLY = True
CONTIG_MAX_DISTANCE_BETWEEN_PAIRS = 1000
# only best
CONTIG_COUNT_MULTI_MAPPED_READS = False
PROKKA_KINGDOM = "Bacteria"

MAXBIN_MAX_ITERATION = 50
MAXBIN_MIN_CONTIG_LENGTH = 200
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

ADAPTERS = "adapters.fa"
RRNA = "silva_rfam_all_rRNAs.fa"
PHIX = "phiX174_virus.fa"

import multiprocessing
import os
import sys
import tempfile

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
