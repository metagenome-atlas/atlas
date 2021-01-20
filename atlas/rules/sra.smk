wildcard_constraints:
    SRR="[S,E,D]RR[0-9]+"

#localrules: prefetch
rule prefetch:
    output:
        sra=touch(temp("SRAreads/{run}_downloaded")),
        # not givins sra file as output allows for continue from the same download
    params:
        outdir='SRAreads'
    log:
        "logs/SRAdownload/{run}.log"
    benchmark:
        "log/benchmarks/logs/SRAdownload/prefetch/{run}.tsv"
    threads:
        1
    resources:
        time= 12, #h
        mem=1,
        internet_connection=1
    conda:
        "%s/sra.yaml" % CONDAENV
    shell:
        " mkdir -p {params.outdir} 2> {log} ;" 
        " cd {params.outdir} 2>> {log};"
        " prefetch "
        " -X 999999999 "
        " {wildcards.run} &>> {log}"


rule extract_run:
    input:
        flag=rules.prefetch.output,
    output:
        "SRAreads/{run}_1.fastq",
        "SRAreads/{run}_2.fastq",
    params:
        outdir='SRAreads',
        sra_file= 'SRAreads/{run}.sra'
    log:
        "logs/SRAdownload/{run}.log"
    threads:
        config['simplejob_threads']
    resources:
        time= lambda wildcards, attempt: attempt * 6, #h
        mem=config['simplejob_mem'],
    conda:
        "%s/sra.yaml" % CONDAENV
    shell:
        "cd {params.outdir} 2>> {log};"
        " fasterq-dump "
        " -e {threads} "
        " &>> {log}"
