import logging
import multiprocessing
import os
import tempfile
import yaml
from collections import OrderedDict
from snakemake.io import load_configfile


ADAPTERS = "adapters.fa"
RRNA = "silva_rfam_all_rRNAs.fa"
PHIX = "phiX174_virus.fa"


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
    """Write the file `config` and complete the sample names and paths for all files in `path`."""
    represent_dict_order = lambda self, data: self.represent_mapping('tag:yaml.org,2002:map', data.items())
    yaml.add_representer(OrderedDict, represent_dict_order)
    path = os.path.realpath(os.path.expanduser(path))

    conf = OrderedDict()
    samples = get_sample_files(path, data_type)
    logging.info("Found %d samples under %s" % (len(samples), path))
    conf["samples"] = samples
    conf["tmpdir"] = tempfile.gettempdir()
    conf["threads"] = multiprocessing.cpu_count() if not threads else threads
    conf["java_mem"] = "16g"
    conf["preprocess_adapters"] = os.path.join(database_dir, ADAPTERS)
    conf["preprocess_adapter_min_k"] = 8
    conf["preprocess_minimum_base_quality"] = 10
    conf["preprocess_allowable_kmer_mismatches"] = 1
    conf["preprocess_reference_kmer_match_length"] = 23
    conf["preprocess_minimum_passing_read_length"] = 51
    conf["preprocess_minimum_base_frequency"] = 0.05

    conf["perform_error_correction"] = "true"

    contamination = OrderedDict()
    contamination["contaminant_references"] = {"rRNA":os.path.join(database_dir, RRNA),
                                               "PhiX":os.path.join(database_dir, PHIX)}
    conf["contaminant_max_indel"] = 20
    conf["contaminant_min_ratio"] = 0.65
    conf["contaminant_kmer_length"] = 12
    conf["contaminant_minimum_hits"] = 1
    conf["contaminant_ambiguous"] = "best"

    conf["normalization_kmer_length"] = 19
    conf["normalization_target_depth"] = 100
    conf["normalization_minimum_kmers"] = 8

    conf["assembler"] = "megahit"
    conf["megahit_memory"] = 0.99
    conf["megahit_min_count"] = 2
    conf["megahit_k_min"] = 21
    conf["megahit_k_max"] = 121
    conf["megahit_k_step"] = 20
    conf["megahit_merge_level"] = "20,0.98"
    conf["megahit_prune_level"] = 2
    conf["megahit_low_local_ratio"] = 0.2
    conf["minimum_contig_length"] = 200
    conf["spades_k"] = "auto"
    conf["minimum_average_coverage"] = 5
    conf["minimum_percent_covered_bases"] = 40
    conf["minimum_mapped_reads"] = 0
    conf["contig_trim_bp"] = 0

    conf["translation_table"] = 11
    conf["minimum_region_overlap"] = 1
    conf["primary_only"] = "false"
    conf["count_multi_mapped_reads"] = "true"
    conf["maximum_counted_map_sites"] = 10
    conf["perform_genome_binning"] = "true"
    conf["maxbin_max_iteration"] = 50
    conf["maxbin_min_contig_length"] = 500
    conf["maxbin_prob_threshold"] = 0.9

    conf["refseq_namemap"] = os.path.join(database_dir, "refseq.db")
    conf["refseq_tree"] = os.path.join(database_dir, "refseq.tree")
    conf["diamond_db"] = os.path.join(database_dir, "refseq.dmnd")
    conf["diamond_run_mode"] = "fast"
    conf["diamond_top_seqs"] = 2
    conf["diamond_e_value"] = 0.000001
    conf["diamond_min_identity"] = 50
    conf["diamond_query_coverage"] = 60
    conf["diamond_gap_open"] = 11
    conf["diamond_gap_extend"] = 1
    conf["diamond_block_size"] = 2
    conf["diamond_index_chunks"] = 1
    conf["summary_method"] = "lca"
    conf["aggregation_method"] = "lca-majority"
    conf["majority_threshold"] = 0.51

    with open(config, "w") as f:
        print(yaml.dump(conf, default_flow_style=False), file=f)
    logging.info("Configuration file written to %s" % config)
