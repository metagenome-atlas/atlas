rule download_SRR_paired:
    output:
        "SRAreads/{SRR}_1.fastq.gz",
        "SRAreads/{SRR}_2.fastq.gz",
        touch(temp('SRAreads/{SRR}_downloaded'))
    params:
        outdir='SRAreads'
    wildcard_constraints:
        SRR="[S,E]RR[0-9]+"
    threads:
        config['simplejob_threads']
    resources:
        time= lambda wildcards, attempt: attempt * 30,
        mem=config['simplejob_mem']
    conda:
        "%s/sra.yaml" % CONDAENV
    shadow:
        "full"
    shell:
        "parallel-fastq-dump --sra-id {wildcards.SRR} --threads {threads} --gzip --split-files --outdir {params.outdir}"

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
