import json
from atlas.utils import split_fasta


rule split:
    input:
        faa = "{sample}/annotation/orfs/{sample}.faa"
    output:
        temp(dynamic("{sample}/annotation/orfs/{sample}_{n}.faa"))
    params:
        chunk_size = config["annotation"].get("chunk_size", 250000)
    run:
        split_fasta(input.faa, chunk_size=params.chunk_size)


rule merge_alignments:
    input:
        dynamic("{sample}/annotation/{reference}/{sample}_intermediate_{n}.aln")
    output:
        "{sample}/annotation/{reference}/{sample}_hits.tsv"
    shell:
        "{SHPFXS} cat {input} | sort -k1,1 -k12,12rn > {output}"


rule parse_blast:
    input:
        "{sample}/annotation/{reference}/{sample}_hits.tsv"
    output:
        "{sample}/annotation/{reference}/{sample}_assignments.tsv"
    params:
        namemap = lambda wc: config["annotation"]["references"][wc.reference]["namemap"],
        treefile = lambda wc: config["annotation"]["references"][wc.reference].get("tree", ""),
        summary_method = lambda wc: config["annotation"]["references"][wc.reference].get("summary_method", "best"),
        aggregation_method = lambda wc: "--aggregation-method %s" % config["annotation"]["references"][wc.reference].get("aggregation_method", "") if "refseq" in wc.reference else "",
        majority_threshold = lambda wc: "--majority-threshold %f" % config["annotation"]["references"][wc.reference].get("majority_threshold", 0.51) if "refseq" in wc.reference else "",
        min_identity = lambda wc: config["annotation"]["references"][wc.reference].get("min_identity", "50"),
        min_bitscore = lambda wc: config["annotation"]["references"][wc.reference].get("min_bitscore", "0"),
        min_length = lambda wc: config["annotation"]["references"][wc.reference].get("min_length", "60"),
        max_evalue = lambda wc: config["annotation"]["references"][wc.reference].get("max_evalue", "0.000001"),
        max_hits = lambda wc: config["annotation"]["references"][wc.reference].get("max_hits", "10"),
        top_fraction = lambda wc: config["annotation"]["references"][wc.reference].get("top_fraction", "0.50")
    shell:
        """{SHPFXS} atlas {wildcards.reference} \
               --summary-method {params.summary_method} {params.aggregation_method} \
               {params.majority_threshold} --min-identity {params.min_identity} \
               --min-bitscore {params.min_bitscore} --min-length {params.min_length} \
               --max-evalue {params.max_evalue} --max-hits {params.max_hits} \
               --top-fraction {params.top_fraction} {input} {params.namemap} {params.treefile} \
               {output}"""


rule merge_blast:
    input:
        ["{sample}/annotation/%s/{sample}_assignments.tsv" % i for i in list(config["annotation"]["references"].keys())]
    output:
        "{sample}/annotation/{sample}_merged_assignments.tsv"
    shell:
        "{SHPFXS} atlas merge-tables {input} {output}"


rule aggregate_counts:
    input:
        merged = "{sample}/annotation/{sample}_merged_assignments.tsv",
        counts = "{sample}/annotation/orfs/{sample}.CDS.txt"
    output:
        ["{sample}/count_tables/{sample}_%s.tsv" % i for i in TABLES]
    params:
        prefix = lambda wc: "{sample}/count_tables/{sample}".format(sample=wc.sample),
        combos = json.dumps(config["summary_counts"])
    shell:
        """{SHPFXS} atlas counts {params.prefix} {input.merged} \
               {input.counts} '{params.combos}'"""
