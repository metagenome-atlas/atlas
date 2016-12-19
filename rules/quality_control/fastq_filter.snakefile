rule quality_filter_reads:
    input:
        lambda wc: config["samples"][wc.sample]["path"]
    output:
        pe = "results/{eid}/{sample}/quality_control/quality_filter/{sample}_pe.fastq.gz",
        se = "results/{eid}/{sample}/quality_control/quality_filter/{sample}_se.fastq.gz",
        stats = "results/{eid}/{sample}/logs/{sample}_quality_filtering_stats.txt"
    params:
        lref = config["preprocessing"]["adapters"],
        rref = config["preprocessing"]["adapters"],
        mink = config["preprocessing"].get("mink", "8"),
        trimq = config["preprocessing"].get("minimum_base_quality", "10"),
        hdist = config["preprocessing"].get("allowable_kmer_mismatches", "1"),
        k = config["preprocessing"].get("reference_kmer_match_length", "31"),
        qtrim = "rl",
        minlength = config["preprocessing"].get("minimum_passing_read_length", "51"),
        minbasefrequency = config["preprocessing"].get("min_base_frequency", 0.05),
        inputs = lambda wc: "in=%s" % config["samples"][wc.sample]["path"][0] if len(config["samples"][wc.sample]["path"]) == 1 else "in=%s in2=%s" % (config["samples"][wc.sample]["path"][0], config["samples"][wc.sample]["path"][1])
    log:
        "results/{eid}/{sample}/logs/{sample}_quality_filter_first_pass.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbduk2.sh {params.inputs} out={output.pe} outs={output.se} \
               rref={params.rref} lref={params.lref} mink={params.mink} qout=33 \
               stats={output.stats} hdist={params.hdist} k={params.k} \
               trimq={params.trimq} qtrim={params.qtrim} threads={threads} \
               minlength={params.minlength} minbasefrequency={params.minbasefrequency} \
               overwrite=true 2> {log}"""
