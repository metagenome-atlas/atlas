wildcard_constraints:
    SRR="[S,E,D]RR[0-9]+"

localrules: prefetch
rule prefetch:
    output:
        sra=temp("SRAreads/{run}.sra"),
    params:
        outdir='SRAreads'
    log:
        "logs/SRAdownload/{run}.log"
    threads:
        1
    shadow:
        'minimal'
    resources:
        time= lambda wildcards, attempt: attempt * 48, #h
        mem=1,
        download=1
    conda:
        "%s/sra.yaml" % CONDAENV
    shell:
        "prefetch -O {params.outdir} "
        "-X 999999999 "
        "{wildcards.SRR} &> {log}"


rule extract_run:
    input:
        sra="SRAreads/{run}.sra",
    output:
        "SRAreads/{run}_1.fastq.gz",
        "SRAreads/{run}_2.fastq.gz",
    params:
        outdir='SRAreads'
    log:
        "logs/SRAdownload/{run}.log"
    threads:
        config['simplejob_threads']
    resources:
        time= lambda wildcards, attempt: attempt * 6, #h
        mem=config['simplejob_mem'],
    conda:
        "%s/sra.yaml" % CONDAENV
    shadow:
        "minimal"
    shell:
        "parallel-fastq-dump "
        "--threads {threads} "
        "--gzip --split-files "
        "--outdir {params.outdir} "
        " -s {input.sra} "
        "&>> {log}"
