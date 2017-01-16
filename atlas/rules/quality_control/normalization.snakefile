rule normalization:
    input:
        "{sample}/quality_control/decontamination/{sample}_pe.fastq.gz"
    output:
        "{sample}/quality_control/{NORMALIZATION}/{sample}_pe.fastq.gz"
    params:
        k = config["preprocessing"]["normalization"].get("k", 21),
        t = config["preprocessing"]["normalization"].get("t", 100),
        minkmers = config["preprocessing"]["normalization"].get("minkmers", 15)
    log:
        "{sample}/logs/{sample}_{NORMALIZATION}.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbnorm.sh in={input} out={output} k={params.k} t={params.t} \
               minkmers={params.minkmers} prefilter=t threads={threads} 2> {log}"""
