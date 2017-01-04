"""
combining samples into a single fastq can take place after intial qc -- all files are interleaved
concatenate into single file; push through the remainder of the workflow
these qc'd files are then used for counts per sample across the assembled contigs
"""
import os
import sys
import tempfile


def get_count_tables(config, key):
    expected_tables = []
    for name, vals in config[key].items():
        if name.lower() == "taxonomy":
            tax_levels = vals.get("levels", ["species"])
            for level in tax_levels:
                level = level.lower()
                tax_name = "taxonomy_%s" % level
                expected_tables.append("%s_%s" % (name, level))
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


def get_temp_dir(config):
    if config.get("tmpdir"):
        return config["tmpdir"]
    else:
        return tempfile.gettempdir()


def coassemblies():
    if "coassemblies" in config["samples"]:
        coassembled_samples = list(config["samples"]["coassemblies"].keys())
        for cs in coassembled_samples:
            for s in config["samples"]["coassemblies"].keys():
                if s not in config["samples"]:
                    print("Coassembly %s includes an undefined sample [%s]" % (cs, s))
                    sys.exit(1)


# rule build_coassembly_fastq:
#     input:
#         lambda wc: ["{sample}/quality_control/decontamination/{sample}_pe.fastq.gz".format(sample=s) for s in config["samples"]["coassemblies"][wc.coassembly]]
#     output:
#         "{coassembly}/reads/{coassembly}_nonorm_pe.fastq.gz"
#     threads:
#         1
#     shell:
#         "cat {input} > {output}"


# rule normalize_coassembly_fastq:
#     "{sample}/quality_control/%s/{sample}_pe.fastq.gz" % NORMALIZATION
#
# "{sample}/quality_control/decontamination/{sample}_pe.fastq.gz"
#
#     coassemblies:
#         test-coa:
#             - shewanella
#             - cytophaga
#             - flavobacterium



if config.get("workflow", "complete") == "complete":

    SHPFXM = config.get("prefix") + str(config.get("threads")) if config.get("prefix") else ""
    SHPFXS = config.get("prefix") + "1" if config.get("prefix") else ""
    SAMPLES = [i for i in config["samples"].keys() if not i == "coassemblies"]
    TABLES = get_count_tables(config, "summary_counts")
    NORMALIZATION = "normalization_k%d_t%d" % (config["preprocessing"]["normalization"].get("k", 21), config["preprocessing"]["normalization"].get("t", 100))
    ASSEMBLER = get_assembler(config) + "_" + NORMALIZATION
    TMPDIR = get_temp_dir(config)

    rule all:
        input:
            expand("{sample}/quality_control/decontamination/{sample}_{decon_dbs}.fastq.gz", sample=SAMPLES, decon_dbs=list(config["preprocessing"]["contamination"]["references"].keys())),
            expand("{sample}/quality_control/decontamination/{sample}_refstats.txt", sample=SAMPLES),
            expand("{sample}/quality_control/%s/{sample}_pe.fastq.gz" % NORMALIZATION, sample=SAMPLES),
            expand("{sample}/quality_control/fastqc/{sample}_pe_fastqc.zip", sample=SAMPLES),
            expand("{sample}/quality_control/fastqc/{sample}_pe_fastqc.html", sample=SAMPLES),
            expand("{sample}/%s/{sample}_contigs.fasta" % ASSEMBLER, sample=SAMPLES),
            expand("{sample}/annotation/orfs/{sample}.faa", sample=SAMPLES),
            expand("{sample}/%s/stats/prefilter_contig_stats.txt" % ASSEMBLER, sample=SAMPLES),
            expand("{sample}/%s/stats/final_contig_stats.txt" % ASSEMBLER, sample=SAMPLES),
            expand("{sample}/annotation/{reference}/{sample}_hits.tsv", reference=list(config["annotation"]["references"].keys()), sample=SAMPLES),
            expand("{sample}/annotation/{reference}/{sample}_assignments.tsv", reference=list(config["annotation"]["references"].keys()), sample=SAMPLES),
            expand("{sample}/annotation/{sample}_merged_assignments.tsv", sample=SAMPLES),
            expand("{sample}/count_tables/{sample}_{table}.tsv", sample=SAMPLES, table=TABLES),
            expand("{sample}/{sample}_readme.html", sample=SAMPLES)

            # coassemblies()

    include: "rules/quality_control/fastq_filter.snakefile"
    include: "rules/quality_control/error_correction.snakefile"
    include: "rules/quality_control/contig_filters.snakefile"
    include: "rules/quality_control/decontamination.snakefile"
    include: "rules/quality_control/normalization.snakefile"
    include: "rules/quality_control/fastqc.snakefile"
    if config["assembly"]["assembler"] == "spades":
        include: "rules/assemblers/spades.snakefile"
    else:
        include: "rules/assemblers/megahit.snakefile"
    include: "rules/annotation/diamond.snakefile"
    include: "rules/annotation/prodigal.snakefile"
    include: "rules/annotation/verse.snakefile"
    include: "rules/annotation/munging.snakefile"
    include: "rules/reports/sample.snakefile"

elif config.get("workflow") == "download":

    FILES = ["silva_rfam_all_rRNAs.fa", "adapters.fa", "phiX174_virus.fa",
             "refseq.db", "refseq.dmnd", "refseq.tree", "cazy.db", "cazy.dmnd",
             "eggnog.db", "eggnog.dmnd", "expazy.db", "expazy.dmnd"]
    OUTDIR = os.path.realpath(config["db_dir"])

    rule all:
        input:
            expand("{dir}/{filename}", dir=OUTDIR, filename=FILES)
    include: "rules/initialization/download.snakefile"

else:
    print("Workflow %s is not a defined workflow." % config.get("workflow", "[no --workflow specified]"), file=sys.stderr)
