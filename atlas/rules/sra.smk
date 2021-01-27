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
        " "
        " vdb-validate {params.outdir}/{wildcards.run}/{wildcards.run}.sra &>> {log} "




rule extract_run:
    input:
        flag=rules.prefetch.output,
    output:
        expand("SRAreads/{{run}}_{fraction}.fastq.gz",
                fraction= ['1','2']
                 )
    params:
        outdir=os.path.abspath('SRAreads'),
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
        " cd {params.outdir} 2>> {log} ;"
        " "
        " vdb-validate {wildcards.run}/{wildcards.run}.sra &>> ../{log} ;"
        " "
        " fasterq-dump "
        " --threads {threads} "
        " --mem {resources.mem}GB "
        " --temp {params.tmpdir}/fasterqdump_tmp/ "
        " --outdir {params.tmpdir}/fasterqdump/ "
        " --log-level debug "
        " --progress "
        " --print-read-nr "
        " {wildcards.run}/{wildcards.run}.sra "
        " &>> ../{log} ; "
        " cd .. ;"
        " "
        " pigz -p{threads} -2 {params.tmpdir}/fasterqdump/{wildcards.run}_?.fastq 2>> {log} ; "
        " mv {params.tmpdir}/fasterqdump/{wildcards.run}_?.fastq.gz "
        "           {params.outdir} 2>> {log} ; "
        " "
        " rm -rf {params.outdir}/{wildcards.run} 2>> {log} ; "
        " rm -f {input} 2>> {log}"

rule download_all_reads:
    input:
        expand("SRAreads/{sample}_{fraction}.fastq.gz",sample=SAMPLES,fraction=['1','2'])
