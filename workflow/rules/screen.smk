
rule generate_sketch:
    input:
        unpack(get_input_fastq),
    output:
        "{sample}/.sketch.gz"
    priority: 100
    log:
        "{sample}/logs/QC/make_sketch.log",
    conda:
        "../enva/required_packages.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_mem"],
        java_mem=int(config["simplejob_mem"] * JAVA_MEM_FRACTION),
    shell:
        "bbsketch.sh "
        "in={input[0]}" # take only one read
        " samplerate=0.5"
        " minkeycount=2 "
        " out= blacklist=nt ssu=f name0=sample1 depth=t overwrite=t"
        "bbsketch.sh in=test_reads/sample2_R1.fastq.gz reads=10M samplerate=0.5 minkeycount=2 out=sample2.sketch blacklist=nt ssu=f name0=sample2 depth=t overwrite=t"


        comparesketch.sh alltoall format=3 out=sketch_comparison.tsv sample?.sketch 

        sendsketch.sh sample2.sketch printdepth2=t level=2 printqfname=f printvolume=t color=f out