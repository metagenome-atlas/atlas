rule split:
    input:
        "results/{eid}/{sample}/annotation/orfs/{sample}.faa"
    output:
        temp(dynamic("results/{eid}/{sample}/annotation/orfs/{sample}_{n}.faa"))
    params:
        chunk_size = config["annotation"].get("chunk_size", "250000")
    shell:
        "{SHPFXS} python scripts/fastx.py split-fasta --chunk-size {params.chunk_size} {input}"


rule merge_alignments:
    input:
        dynamic("results/{eid}/{sample}/annotation/{reference}/{sample}_intermediate_{n}.aln")
    output:
        "results/{eid}/{sample}/annotation/{reference}/{sample}_hits.tsv"
    shell:
        "{SHPFXS} cat {input} | sort -k1,1 -k12,12rn > {output}"


rule parse_blast:
    input:
        "results/{eid}/{sample}/annotation/{reference}/{sample}_hits.tsv"
    output:
        "results/{eid}/{sample}/annotation/{reference}/{sample}_assignments.tsv"
    params:
        # subcommand = lambda wc: "refseq" if "refseq" in wc.reference else "eggnog",
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
        """{SHPFXS} python scripts/blast2assignment.py {wildcards.reference} \
               --summary-method {params.summary_method} {params.aggregation_method} \
               {params.majority_threshold} --min-identity {params.min_identity} \
               --min-bitscore {params.min_bitscore} --min-length {params.min_length} \
               --max-evalue {params.max_evalue} --max-hits {params.max_hits} \
               --top-fraction {params.top_fraction} {input} {params.namemap} {params.treefile} \
               {output}"""


rule merge_blast:
    input:
        ["results/{eid}/{sample}/annotation/%s/{sample}_assignments.tsv" % i for i in list(config["annotation"]["references"].keys())]
    output:
        "results/{eid}/{sample}/annotation/{sample}_merged_assignments.tsv"
    shell:
        "{SHPFXS} python scripts/blast2assignment.py merge-tables {input} {output}"


rule aggregate_counts:
    input:
        merged = "results/{eid}/{sample}/annotation/{sample}_merged_assignments.tsv",
        counts = "results/{eid}/{sample}/annotation/orfs/{sample}.CDS.txt"
    output:
        ["results/{eid}/{sample}/count_tables/{sample}_%s.tsv" % i for i in TABLES]
    params:
        prefix = lambda wc: "results/%s/%s/count_tables/%s" % (wc.eid, wc.sample, wc.sample),
        combos = json.dumps(config["summary_counts"])
    shell:
        """{SHPFXS} python scripts/blast2assignment.py counts {params.prefix} {input.merged} \
               {input.counts} '{params.combos}'"""
