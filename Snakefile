import json
import os
import sys
from glob import glob
from snakemake.utils import report
from subprocess import check_output


def read_count(fastq, out=None):
    total = 0
    count_file = fastq + ".count" if out is None else out
    if os.path.exists(fastq) and os.path.getsize(fastq) > 100:
        if not os.path.exists(count_file):
            if fastq.endswith("gz"):
                check_output("gunzip -cd %s | awk '{n++}END{print n/4}' > %s" % (fastq, fastq + '.count'), shell=True)
            else:
                check_output("awk '{n++}END{print n/4}' %s > %s" % (fastq, fastq + '.count'), shell=True)

        with open(count_file) as fh:
            for line in fh:
                total = int(line.strip())
                break
    return total


def get_samples(config, coverage_cutoff=1000):
    samples = set()
    omitted = set()
    for sample_id, vals in config["samples"].items():
        count = 0
        for fastq in vals["path"]:
            if len(vals["path"]) == 1:
                count = read_count(fastq) / 2
            else:
                count = read_count(fastq)
            # only need a count for one of the fastqs
            break
        if count > coverage_cutoff:
            samples.add(sample_id)
        else:
            print("Omitting %s due to low starting read count" % sample_id, file=sys.stderr)
            omitted.add(sample_id)
    return samples, omitted


def get_count_tables(config, key):
    expected_tables = []
    for name, vals in config[key].items():
        if name.lower() == "taxonomy":
            tax_levels = vals.get("levels", ["species"])
            for level in tax_levels:
                level = level.lower()
                tax_name = "taxonomy_%s" % level
                for subname, subvals in vals.items():
                    if subname.lower() == "levels": continue
                    expected_tables.append("%s_%s" % (subname, tax_name))
        else:
            expected_tables.append(name)
    return expected_tables


def get_assembler(config):
    if config["assembly"]["assembler"] == "spades":
        return "spades_{k}".format(k=config["assembly"].get("spades_k", "auto").replace(",", "_"))
    else:
        k_min = config["assembly"].get("kmer_min", 21)
        k_max = config["assembly"].get("kmer_max", 121)
        k_step = config["assembly"].get("kmer_step", 20)
        return "megahit_{min}_{max}_{step}".format(min=k_min, max=k_max, step=k_step)


# shell prefixes for multi-threaded and single-threads tasks
SHPFXM = config.get("prefix") + str(config.get("threads")) if config.get("prefix") else ""
SHPFXS = config.get("prefix") + "1" if config.get("prefix") else ""
SAMPLES, OMITTED = get_samples(config)
TABLES = get_count_tables(config, "summary_counts")
NORMALIZATION = "normalization_k%d_t%d" % (config["preprocessing"]["normalization"].get("k", 31), config["preprocessing"]["normalization"].get("t", 100))
ASSEMBLER = get_assembler(config)


rule all:
    input:
        expand("results/{eid}/{sample}/quality_control/decontamination/{sample}_{decon_dbs}.fastq.gz", eid=config["experiment"], sample=SAMPLES, decon_dbs=list(config["preprocessing"]["contamination"]["references"].keys())),
        expand("results/{eid}/{sample}/quality_control/decontamination/{sample}_refstats.txt", eid=config["experiment"], sample=SAMPLES),
        expand("results/{eid}/{sample}/quality_control/%s/{sample}_pe.fastq.gz" % NORMALIZATION, eid=config["experiment"], sample=SAMPLES),
        expand("results/{eid}/{sample}/quality_control/fastqc/{sample}_final_fastqc.zip", eid=config["experiment"], sample=SAMPLES),
        expand("results/{eid}/{sample}/quality_control/fastqc/{sample}_final_fastqc.html", eid=config["experiment"], sample=SAMPLES),
        expand("results/{eid}/{sample}/%s/{sample}_contigs.fasta" % ASSEMBLER, eid=config["experiment"], sample=SAMPLES),
        expand("results/{eid}/{sample}/annotation/orfs/{sample}.faa", eid=config["experiment"], sample=SAMPLES),
        expand("results/{eid}/{sample}/%s/stats/prefilter_contig_stats.txt" % ASSEMBLER, eid=config["experiment"], sample=SAMPLES),
        expand("results/{eid}/{sample}/%s/stats/final_contig_stats.txt" % ASSEMBLER, eid=config["experiment"], sample=SAMPLES),
        expand("results/{eid}/{sample}/annotation/{reference}/{sample}_hits.tsv", eid=config["experiment"], reference=list(config["annotation"]["references"].keys()), sample=SAMPLES),
        expand("results/{eid}/{sample}/annotation/{reference}/{sample}_assignments.tsv", eid=config["experiment"], reference=list(config["annotation"]["references"].keys()), sample=SAMPLES),
        expand("results/{eid}/{sample}/annotation/{sample}_merged_assignments.tsv", eid=config["experiment"], sample=SAMPLES),
        expand("results/{eid}/{sample}/count_tables/{sample}_{table}.tsv", eid=config["experiment"], sample=SAMPLES, table=TABLES),
        expand("results/{eid}/{sample}/{sample}_readme.html", eid=config["experiment"], sample=SAMPLES)


if config["assembly"]["assembler"] == "spades":
    include: "rules/assemblers/spades.rules"
else:
    include: "rules/assemblers/megahit.rules"

include: "rules/quality_control/fastq_filter.rules"
include: "rules/quality_control/error_correction.rules"
include: "rules/quality_control/contig_filters.rules"
include: "rules/quality_control/decontamination.rules"
include: "rules/quality_control/normalization.rules"
include: "rules/quality_control/fastqc.rules"
include: "rules/annotation/diamond.rules"
include: "rules/annotation/prodigal.rules"
include: "rules/annotation/munging.rules"
include: "rules/reports/sample.rules"
