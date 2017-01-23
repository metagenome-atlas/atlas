import logging
import multiprocessing
import os
import tempfile
import yaml
from collections import OrderedDict


ADAPTERS = "adapters.fa"
RRNA = "silva_rfam_all_rRNAs.fa"
PHIX = "phiX174_virus.fa"

CAZY = "cazy"
COG = "cog"
EGGNOG = "eggnog"
EXPAZY = "expazy"
REFSEQ = "refseq"


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
                samples[sample_id] = {'path': fastq_paths, 'type': data_type}
    return samples


def make_config(config, path, data_type, database_dir, threads, assembler):
    """Write the file `config` and complete the sample names and paths for all files in `path`."""
    represent_dict_order = lambda self, data:  self.represent_mapping('tag:yaml.org,2002:map', data.items())
    yaml.add_representer(OrderedDict, represent_dict_order)
    path = os.path.realpath(path)

    conf = OrderedDict()
    samples = get_sample_files(path, data_type)
    logging.info("Found %d samples under %s" % (len(samples), path))
    conf["samples"] = samples
    conf["tmpdir"] = tempfile.gettempdir()
    conf["database_directory"] = database_dir
    conf["threads"] = multiprocessing.cpu_count() if not threads else threads
    # conf["prefix"]

    preprocessing = OrderedDict()
    preprocessing["adapters"] = os.path.join(database_dir, ADAPTERS)

    preprocessing["mink"] = 8
    preprocessing["minimum_base_quality"] = 10
    preprocessing["allowable_kmer_mismatches"] = 1
    preprocessing["reference_kmer_match_length"] = 31
    # preprocessing["qtrim"] = "rl"
    preprocessing["minimum_passing_read_length"] = 51
    preprocessing["min_base_frequency"] = 0.05

    contamination = OrderedDict()
    contamination["references"] = {"rRNA":os.path.join(database_dir, RRNA), "PhiX":os.path.join(database_dir, PHIX)}
    contamination["maxindel"] = 20
    contamination["minratio"] = 0.65
    contamination["minhits"] = 1
    contamination["ambiguous"] = "best"
    contamination["k"] = 15

    preprocessing["contamination"] = contamination
    preprocessing["normalization"] = {"k":21, "t":100, "minkmers":8}
    conf["preprocessing"] = preprocessing

    assembly = OrderedDict()
    assembly["assembler"] = assembler
    # assembly["memory"] = 0.99
    assembly["minimum_count"] = 2
    assembly["kmer_min"] = 21
    assembly["kmer_max"] = 121
    assembly["kmer_step"] = 20
    # assembly["merge_level"] = "20,0.98"
    # assembly["prune_level"] = 2
    # assembly["low_local_ratio"] = 0.2
    assembly["minimum_contig_length"] = 200
    assembly["spades_k"] = "auto"
    assembly["minc"] = 5
    assembly["minp"] = 40
    assembly["minr"] = 0
    assembly["minl"] = 250
    assembly["trim"] = 0
    conf["assembly"] = assembly

    annotation = OrderedDict()
    # annotation["translation_table"]
    annotation["minimum_overlap"] = 20
    eggnog = OrderedDict()
    eggnog["namemap"] = os.path.join(database_dir, "%s.db" % EGGNOG)
    eggnog["dmnd"] = os.path.join(database_dir, "%s.dmnd" % EGGNOG)
    eggnog["run_mode"] = "fast"
    eggnog["top_seqs"] = 5
    eggnog["summary_method"] = "best"
    annotation["eggnog"] = eggnog

    refseq = OrderedDict()
    refseq["namemap"] = os.path.join(database_dir, "%s.db" % REFSEQ)
    refseq["tree"] = os.path.join(database_dir, "%s.tree" % REFSEQ)
    refseq["dmnd"] = os.path.join(database_dir, "%s.dmnd" % REFSEQ)
    refseq["run_mode"] = "fast"
    refseq["top_seqs"] = 5
    refseq["summary_method"] = "best"
    refseq["aggregation_method"] = "lca-majority"
    refseq["majority_threshold"] = 0.51
    annotation["refseq"] = refseq

    expazy = OrderedDict()
    expazy["namemap"] = os.path.join(database_dir, "%s.db" % EXPAZY)
    expazy["dmnd"] = os.path.join(database_dir, "%s.dmnd" % EXPAZY)
    expazy["run_mode"] = "fast"
    expazy["top_seqs"] = 2
    expazy["summary_method"] = "majority"
    expazy["index_chunks"] = 1
    annotation["expazy"] = expazy

    cazy = OrderedDict()
    cazy["namemap"] = os.path.join(database_dir, "%s.db" % CAZY)
    cazy["dmnd"] = os.path.join(database_dir, "%s.dmnd" % CAZY)
    cazy["run_mode"] = "fast"
    cazy["top_seqs"] = 2
    cazy["summary_method"] = "majority"
    cazy["index_chunks"] = 1
    annotation["cazy"] = cazy

    cog = OrderedDict()
    cog["namemap"] = os.path.join(database_dir, "%s.db", COG)
    cog["dmnd"] = os.path.join(database_dir, "%s.dmnd" % COG)
    cog["run_mode"] = "fast"
    cog["top_seqs"] = 2
    cog["summary_method"] = "majority"
    cog["index_chunks"] = 1
    annotation["cog"] = cog

    conf["annotation"] = annotation

    summary_counts = OrderedDict()

    summary_counts["taxonomy"] = {"levels":["phylum", "class", "order", "species"],
                                  "CAZy_Family":["cazy_family"],
                                  "ExPAZy":["expazy_name", "expazy_ec"],
                                  "RefSeq":["refseq_product"],
                                  "COG":["cog_id", "cog_functional_class", "cog_annotation"]}
    summary_counts["KO"] = ["ko_id", "ko_gene_symbol", "ko_product", "ko_ec"]
    summary_counts["RefSeq"] = ["refseq_product"]
    summary_counts["COG"] = ["cog_id", "cog_functional_class", "cog_annotation"]
    summary_counts["ExPAZy"] = ["expazy_name", "expazy_ec"]
    summary_counts["CAZy"] = ["cazy_family", "cazy_class"]

    conf["summary_counts"] = summary_counts

    with open(config, "w") as f:
        print(yaml.dump(conf, default_flow_style=False), file=f)
    logging.info("Configuration file written to %s" % config)
