



rule identify:
    input:
        dir=genome_dir,
        flag= rules.download_gtdb.output
    output:
        directory("genomes/taxonomy/identify")
    threads:
        config['threads']
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk_identify.txt",
        "genomes/taxonomy/gtdbtk.log"
    params:
        outdir= "genomes/taxonomy",
        extension="fasta",
    shell:
        "GTDBTK_DATA_PATH={GTDBTK_DATA_PATH} ; "
        "gtdbtk identify --genome_dir {input.dir} --out_dir {params.outdir} "
        "--extension {params.extension} "
        "--cpus {threads} &> {log[0]}"

rule align:
    input:
        "genomes/taxonomy/identify"
    output:
        "genomes/taxonomy/gtdbtk.bac120.user_msa.fasta"
    threads:
        config['threads']
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk_align.txt",
        "genomes/taxonomy/gtdbtk.log"
    params:
        outdir= "genomes/taxonomy"
    shell:
        "GTDBTK_DATA_PATH={GTDBTK_DATA_PATH} ;  "
        "gtdbtk align --identify_dir {params.outdir} --out_dir {params.outdir} "
        "--cpus {threads} &> {log[0]}"


rule classify:
    input:
        rules.align.output,
        genome_dir=genome_dir,
    output:
        "genomes/taxonomy/gtdbtk.bac120.summary.tsv",
    threads:
        8 #pplacer needs much memory for not many threads
    resources:
        mem=conf['large_mem']
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk_classify.txt",
        "genomes/taxonomy/gtdbtk.log"
    params:
        outdir= "genomes/taxonomy",
        extension="fasta",
    shell:
        "GTDBTK_DATA_PATH={GTDBTK_DATA_PATH} ; "
        "gtdbtk classify --genome_dir {input.genome_dir} --align_dir {params.outdir} "
        "--out_dir {params.outdir} "
        "--extension {params.extension} "
        "--cpus {threads} &> {log[0]}"


rule infer:
    input:
        "genomes/taxonomy/gtdbtk.bac120.user_msa.fasta"
    output:
        "genomes/taxonomy/gtdbtk.unrooted.tree"
    threads:
        config['threads']
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk_infer.txt",
        "genomes/taxonomy/gtdbtk.log"
    params:
        outdir="genomes/taxonomy"
    shell:
        "GTDBTK_DATA_PATH={GTDBTK_DATA_PATH} ;  "
        "gtdbtk infer --msa_file {input} --out_dir {params.outdir} "
        "--cpus {threads} &> {log[0]}"
