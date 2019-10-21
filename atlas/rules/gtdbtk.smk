

gtdb_dir="genomes/taxonomy/gtdb"

rule identify:
    input:
        dir=genome_dir,
        flag= rules.download_gtdb.output
    output:
        directory(f"{gtdb_dir}/identify")
    threads:
        config['threads']
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk/identify.txt",
        f"{gtdb_dir}/gtdbtk.log"
    params:
        outdir= gtdb_dir,
        extension="fasta",
    shell:
        "GTDBTK_DATA_PATH={GTDBTK_DATA_PATH} ; "
        "gtdbtk identify --genome_dir {input.dir} --out_dir {params.outdir} "
        "--extension {params.extension} "
        "--cpus {threads} &> {log[0]}"

checkpoint align:
    input:
        f"{gtdb_dir}/identify"
    output:
        directory(f"{gtdb_dir}/align")
    threads:
        config['threads']
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk/align.txt",
        f"{gtdb_dir}/gtdbtk.log"
    params:
        outdir= gtdb_dir
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
        config['threads'] #pplacer needs much memory for not many threads
    resources:
        mem=config['large_mem']
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk/classify.txt",
        f"{gtdb_dir}/gtdbtk.log"
    params:
        outdir= gtdb_dir,
        extension="fasta",
    shell:
        "GTDBTK_DATA_PATH={GTDBTK_DATA_PATH} ; "
        "gtdbtk classify --genome_dir {input.genome_dir} --align_dir {params.outdir} "
        "--out_dir {params.outdir} "
        "--extension {params.extension} "
        "--cpus {threads} &> {log[0]}"


rule infer:
    input:
        f"{gtdb_dir}/gtdbtk.bac120.user_msa.fasta"
    output:
        f"{gtdb_dir}/gtdbtk.unrooted.tree"
    threads:
        config['threads']
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk/infer.txt",
        f"{gtdb_dir}/gtdbtk.log"
    params:
        outdir=gtdb_dir
    shell:
        "GTDBTK_DATA_PATH={GTDBTK_DATA_PATH} ;  "
        "gtdbtk infer --msa_file {input} --out_dir {params.outdir} "
        "--cpus {threads} &> {log[0]}"
