


rule bbsplit_index:
    input:
        genome_dir,
    output:
        directory("ref/genome/5"),
        directory("ref/index/5"),
    params:
        path=lambda wc, output: os.path.dirname(output[0]),
    threads: config["threads"]
    resources:
        mem_mb=config["mem"] * 1000,
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
        time_min=60 * config["runtime"]["long"],
    log:
        "logs/genomes/mapping/build_bbmap_index.log",
    benchmark:
        "logs/benchmark/genomes/bbsplit_index.log"
    conda:
        "../envs/required_packages.yaml"
    shell:
        " bbsplit.sh ref={input[0]} "
        " build=5 "
        " threads={threads} "
        " log={log} "
        " -Xmx{resources.java_mem}G "
        " 2> {log} "


rule bbsplit:
    input:
        reads=get_quality_controlled_reads,
        refdir=rules.bbsplit_index.output,
    output:
        scafstats=temp("genomes/alignments/scafstats/{sample}.tsv.gz"),
        refstats=temp("genomes/alignments/refstats/{sample}.tsv.gz"),
        covstats="genomes/alignments/coverage/{sample}.tsv.gz",
        bincov=temp("genomes/alignments/coverage_binned/{sample}.tsv.gz"),
    params:
        reads=lambda wc, input: io_params_for_tadpole(input.reads, "in"),
        ambiguous="all",
        ambiguous2="best",
        minid=config["contig_min_id"],
        # unmapped=lambda wc, output: io_params_for_tadpole(output.unmapped, "outu"),
    shadow:
        "shallow"
    log:
        "logs/genomes/bbsplit/{sample}.log",
    benchmark:
        "logs/benchmark/genomes/bbsplit/{sample}.txt"
    threads: config["threads"]
    conda:
        "../envs/required_packages.yaml"
    resources:
        mem_mb=1000 * config["mem"],
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
        time_min=config["runtime"]["default"] * 60,
        time=config["runtime"]["default"],
    shell:
        " bbsplit.sh "
        " {params.reads} "
        " build=5 "
        " threads={threads} "
        " ambiguous={params.ambiguous} "
        " ambiguous2={params.ambiguous2} "
        " minid={params.minid} "
        " scafstats={output.scafstats} "
        " refstats={output.refstats} "
        " covstats={output.covstats} "
        " bincov={output.bincov} "
        " 2> {log} "


localrules:
    count_mapped_reads,


rule count_mapped_reads:
    input:
        logfile=expand("logs/genomes/bbsplit/{sample}.log", sample=SAMPLES),
        out_check=expand(rules.bbsplit.output.covstats, sample=SAMPLES),
    output:
        "genomes/counts/mapping_rate.tsv",
    run:
        import pandas as pd
        from utils.parsers_bbmap import parse_bbmap_log_file

        D = pd.DataFrame(index=SAMPLES, columns=["reads_used", "reads_mapped"])
        for i, sample in enumerate(SAMPLES):
            D.loc[sample] = parse_bbmap_log_file(input[i])

        D["mapping_rate"] = D.iloc[:, 1] / D.iloc[:, 0]
        D.to_csv(output[0], sep="\t")


rule combine_bined_coverages_MAGs:
    input:
        binned_coverage_files=expand(
            "genomes/alignments/coverage_binned/{sample}.tsv.gz", sample=SAMPLES
        ),
    params:
        samples=SAMPLES,
    output:
        binned_cov="genomes/counts/binned_coverage.parquet",
        median_abund="genomes/counts/median_coverage_genomes.parquet",
    log:
        "logs/genomes/counts/combine_binned_coverages_MAGs.log",
    threads: 1
    resources:
        mem_mb=1000 * config["simplejob_mem"],
        time_min=config["runtime"]["simplejob"] * 60,
    script:
        "../scripts/combine_coverage_MAGs.py"


rule merge_counts:
    input:
        expand("genomes/alignments/{{scope}}stats/{sample}.tsv.gz", sample=SAMPLES),
    output:
        counts="genomes/counts/counts_{scope}.parquet",
    threads: 1
    run:
        import pandas as pd
        import gc

        samples = SAMPLES

        Reads = {}
        for i in range(len(input)):
            df = pd.read_csv(input[i], index_col=0, sep="\t")
            Reads[samples[i]] = pd.to_numeric(df.unambiguousReads)
            del df
            gc.collect()


        Reads = pd.concat(Reads, axis=1, sort=False).fillna(0).T
        gc.collect()

        Reads.reset_index().to_parquet(output.counts)
