
rule generate_sketch:
    input:
        unpack(get_input_fastq),
    output:
        "Intermediate/screen/sketches/{sample}.sketch.gz"
    log:
        "logs/screen/make_sketch/{sample}.log",
    conda:
        "../envs/required_packages.yaml"
    threads: 1
    resources:
        mem=config["simplejob_mem"],
        java_mem=int(config["simplejob_mem"] * JAVA_MEM_FRACTION),
    shell:
        "bbsketch.sh "
        "in={input[0]}" # take only one read
        " samplerate=0.5"
        " minkeycount=2 "
        " out={output} "
        " blacklist=nt ssu=f name0={wildcards.sample} depth=t overwrite=t "
        " -Xmx{resources.java_mem}g "
        " &> {log}"


rule compare_sketch:
    input:
        expand( rules.generate_sketch.output, sample =SAMPLES ),
    output:
        "QC/screen/sketch_comparison.tsv.gz"
    priority: 100
    log:
        "logs/screen/compare_sketch.log",
    conda:
        "../envs/required_packages.yaml"
    threads: 1
    resources:
        mem=config["mem"],
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
    shell:
        "comparesketch.sh alltoall "
        " format=3 out={output} "
        " {input} "
        " -Xmx{resources.java_mem}g "
        " &> {log}"



#        sendsketch.sh sample2.sketch printdepth2=t level=2 printqfname=f printvolume=t color=f out