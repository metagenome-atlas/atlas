import logging
import multiprocessing
import os
import tempfile
import yaml
from collections import OrderedDict
from snakemake.io import load_configfile

from .default_values import *



def get_sample_files(path, data_type):
    samples = OrderedDict()
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
    """Write the file `config` and complete the sample names and paths for all
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

    represent_dict_order = lambda self, data: self.represent_mapping('tag:yaml.org,2002:map', data.items())
    yaml.add_representer(OrderedDict, represent_dict_order)
    path = os.path.realpath(os.path.expanduser(path))
    database_dir = os.path.realpath(os.path.expanduser(database_dir))

    conf = OrderedDict()
    samples = get_sample_files(path, data_type)
    logging.info("Found %d samples under %s" % (len(samples), path))
    conf["samples"] = samples
    conf["tmpdir"] = tempfile.gettempdir()
    conf["threads"] = multiprocessing.cpu_count() if not threads else threads
    conf["java_mem"] = JAVA_MEM
    conf["preprocess_adapters"] = os.path.join(database_dir, ADAPTERS)
    conf["preprocess_adapter_min_k"] = PREPROCESS_ADAPTER_MIN_K
    conf["preprocess_minimum_base_quality"] = PREPROCESS_MINIMUM_BASE_QUALITY
    conf["preprocess_allowable_kmer_mismatches"] = PREPROCESS_ALLOWABLE_KMER_MISMATCHES
    conf["preprocess_reference_kmer_match_length"] = PREPROCESS_REFERENCE_KMER_MATCH_LENGTH
    conf["preprocess_minimum_passing_read_length"] = PREPROCESS_MINIMUM_PASSING_READ_LENGTH
    conf["preprocess_minimum_base_frequency"] = PREPROCESS_MINIMUM_BASE_FREQUENCY

    conf["deduplicate"] = True
    conf["error_correction_overlapping_pairs"] = True
    conf["assembly_preprocessing_steps"]=['normalized','errorcorr','merged']

    conf["contaminant_references"] = {"rRNA":os.path.join(database_dir, RRNA),
                                      "PhiX":os.path.join(database_dir, PHIX)}
    conf["contaminant_max_indel"] = CONTAMINANT_MAX_INDEL
    conf["contaminant_min_ratio"] = CONTAMINANT_MIN_RATIO
    conf["contaminant_kmer_length"] = CONTAMINANT_KMER_LENGTH
    conf["contaminant_minimum_hits"] = CONTAMINANT_MINIMUM_HITS
    conf["contaminant_ambiguous"] = CONTAMINANT_AMBIGUOUS

    conf["normalization_kmer_length"] = NORMALIZATION_KMER_LENGTH
    conf["normalization_target_depth"] = NORMALIZATION_TARGET_DEPTH
    conf["normalization_minimum_kmers"] = NORMALIZATION_MINIMUM_KMERS

    conf["merging_k"]= MERGING_K
    conf["merging_extend2"] = MERGING_EXTEND2
    conf["merging_flags"]  = MERGING_FLAGS

    conf["assembler"] = assembler
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
    conf["spades_k"] = SPADES_K
    conf["minimum_average_coverage"] = MINIMUM_AVERAGE_COVERAGE
    conf["minimum_percent_covered_bases"] = MINIMUM_PERCENT_COVERED_BASES
    conf["minimum_mapped_reads"] = MINIMUM_MAPPED_READS
    conf["contig_trim_bp"] = CONTIG_TRIM_BP

    conf["translation_table"] = 11
    conf["minimum_region_overlap"] = MINIMUM_REGION_OVERLAP
    conf["primary_only"] = False
    conf["count_multi_mapped_reads"] = True
    conf["maximum_counted_map_sites"] = MAXIMUM_COUNTED_MAP_SITES
    conf["perform_genome_binning"] = True
    conf["maxbin_max_iteration"] = MAXBIN_MAX_ITERATION
    conf["maxbin_min_contig_length"] = MAXBIN_MIN_CONTIG_LENGTH
    conf["maxbin_prob_threshold"] = MAXBIN_PROB_THRESHOLD

    conf["refseq_namemap"] = os.path.join(database_dir, "refseq.db")
    conf["refseq_tree"] = os.path.join(database_dir, "refseq.tree")
    conf["diamond_db"] = os.path.join(database_dir, "refseq.dmnd")
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

    with open(config, "w") as f:
        print(yaml.dump(conf, default_flow_style=False), file=f)
    logging.info("Configuration file written to %s" % config)
