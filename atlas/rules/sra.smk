wildcard_constraints:
    SRR="[S,E,D]RR[0-9]+"

localrules: prefetch
rule prefetch:
    output:
        sra=temp(touch("SRAreads/{run}_downloaded")),
        # not givins sra file as output allows for continue from the same download
    params:
        outdir= 'SRAreads' #lambda wc,output: os.path.dirname(output[0])
    log:
        "log/SRAdownload/{run}.log"
    benchmark:
        "log/benchmarks/SRAdownload/prefetch/{run}.tsv"
    threads:
        1
    resources:
        mem=1,
        time= int(config["runtime"]["simple_job"]),
        internet_connection=1
    conda:
        "%s/sra.yaml" % CONDAENV
    shell:
        " mkdir -p {params.outdir} 2> {log}; "
        " prefetch "
        " --output-directory {params.outdir} "
        " -X 999999999 "
        " --progress "
        " --log-level info "
        " {wildcards.run} &>> {log}"

ruleorder: gzip > extract_run
rule gzip:
    input:
        "SRAreads/{run}_{direction}.fastq",
    output:
        "SRAreads/{run}_{direction}.fastq.gz",
    threads:
        config['simplejob_threads']
    resources:
        time= int(config["runtime"]["simple_job"]),
        mem=1 #default 100Mb
    conda:
        "%s/sra.yaml" % CONDAENV
    shell:
        " pigz -p{threads} -2 {input}"


rule extract_run:
    input:
        flag=rules.prefetch.output,
    output:
        "SRAreads/{run}_1.fastq.gz",
        "SRAreads/{run}_2.fastq.gz",
    params:
        outdir='SRAreads',
        tmpdir= TMPDIR
    log:
        "log/SRAdownload/{run}.log"
    benchmark:
        "log/benchmarks/SRAdownload/fasterqdump/{run}.tsv"
    threads:
        config['simplejob_threads']
    resources:
        time= int(config["runtime"]["simple_job"]),
        mem=1 #default 100Mb
    conda:
        "%s/sra.yaml" % CONDAENV
    shell:
        " cd {params.outdir} 2>> {log};"
        " "
        " fasterq-dump "
        " --threads {threads} "
        " --mem{resources.mem}GB "
        " --temp {params.tmpdir} "
        " --outdir {params.tmpdir} "
        " --log-level info "
        " --progress "
        " --print-read-nr "
        " {wildcards.run} "
        " &>> ../{log} ; "
        " "
        " pigz -p{threads} -2 {params.tmpdir}/{wildcards.run}_1.fastq "
        " mv {params.tmpdir}/{wildcards.run}_1.fastq {output[0]}"
        " "
        " pigz -p{threads} -2 {params.tmpdir}/{wildcards.run}_2.fastq "
        " mv {params.tmpdir}/{wildcards.run}_2.fastq {output[1]}"
        " "
        " rm -rf {wildcards.run} 2> ../{log} "

rule download_all_reads:
    input:
        expand("SRAreads/{sample}_{fraction}.fastq.gz",sample=SAMPLES,fraction=['1','2'])
