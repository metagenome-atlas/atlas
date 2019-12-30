localrules: prefetch
rule prefetch:
    output:
        sra=temp("SRAreads/{SRR}.sra"),
    params:
        outdir='SRAreads'
    wildcard_constraints:
        SRR="[S,E]RR[0-9]+"
    log:
        "logs/SRAdownload/{SRR}.log"
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
        '--ascp-path "/home/kiesers/.aspera/connect/bin/ascp|/home/kiesers/.aspera/connect/etc/asperaweb_id_dsa.openssh" '
        "{wildcards.SRR} &> {log}"


rule download_SRR_paired:
    input:
        sra="SRAreads/{SRR}.sra",
    output:
        "SRAreads/{SRR}_1.fastq.gz",
        "SRAreads/{SRR}_2.fastq.gz",
    params:
        outdir='SRAreads'
    wildcard_constraints:
        SRR="[S,E]RR[0-9]+"
    log:
        "logs/SRAdownload/{SRR}.log"
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

localrules: rename_SRR
rule rename_SRR:
    input:
        rules.download_SRR_paired.output
    output:
        temp("SRAreads/{SRR}_R1.fastq.gz"),
        temp("SRAreads/{SRR}_R2.fastq.gz")
    threads:
        1
    shell:
        "mv {input[0]} {output[0]}; "
        "mv {input[1]} {output[1]}"
