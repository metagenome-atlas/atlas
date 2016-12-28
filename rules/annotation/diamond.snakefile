rule build_dmnd_database:
    input:
        lambda wc: config["annotation"]["references"][wc.reference]["fasta"]
    output:
        "%s/{reference}.dmnd" % config.get("database_directory", "databases")
    threads:
        config.get("threads", 1)
    shell:
        "{SHPFXM} diamond makedb --no-auto-append --threads {threads} --in {input} --db {output}"


rule diamond_alignments:
    input:
        fasta = "{sample}/annotation/orfs/{sample}_{n}.faa",
        db = "%s/{reference}.dmnd" % config.get("database_directory", "databases")
    output:
        temp("{sample}/annotation/{reference}/{sample}_intermediate_{n}.aln")
    params:
        tmpdir = "--tmpdir %s" % TMPDIR if TMPDIR else "",
        top_seqs = lambda wc: config["annotation"]["references"][wc.reference].get("top_seqs", "5"),
        e_value = lambda wc: config["annotation"]["references"][wc.reference].get("e_value", "0.000001"),
        min_identity = lambda wc: config["annotation"]["references"][wc.reference].get("min_identity", "50"),
        query_cover = lambda wc: config["annotation"]["references"][wc.reference].get("query_coverage", "60"),
        gap_open = lambda wc: config["annotation"]["references"][wc.reference].get("gap_open", "11"),
        gap_extend = lambda wc: config["annotation"]["references"][wc.reference].get("gap_extend", "1"),
        block_size = lambda wc: config["annotation"]["references"][wc.reference].get("block_size", "2"),
        index_chunks = lambda wc: config["annotation"]["references"][wc.reference].get("index_chunks", "4"),
        run_mode = lambda wc: "" if config["annotation"]["references"][wc.reference].get("run_mode", "fast") == "fast" else "--more-sensitive"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} diamond blastp --threads {threads} --outfmt 6 --out {output} \
               --query {input.fasta} --db {input.db} --top {params.top_seqs} \
               --evalue {params.e_value} --id {params.min_identity} \
               --query-cover {params.query_cover} {params.run_mode} --gapopen {params.gap_open} \
               --gapextend {params.gap_extend} {params.tmpdir} --block-size {params.block_size} \
               --index-chunks {params.index_chunks}"""
