

GTDBTK_DATA_PATH="/home/kiesers/scratch/Atlas/databases/GTDB-TK"



rule all:
    input:
        "taxonomy/gtdbtk.bac120.summary.tsv",
        "taxonomy/gtdbtk.unrooted.tree"

rule identify:
    input:
        genome_dir
    output:
        directory("taxonomy/identify")
    threads:
        config['threads']
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk_identify.txt",
        "taxonomy/gtdbtk.log"
    params:
        outdir= "taxonomy",
        extension="fasta",
    shell:
        "GTDBTK_DATA_PATH={GTDBTK_DATA_PATH} ; "
        "gtdbtk identify --genome_dir {input} --out_dir {params.outdir} "
        "--extension {params.extension} "
        "--cpus {threads} &> {log[0]}"

rule align:
    input:
        "taxonomy/identify"
    output:
        "taxonomy/gtdbtk.bac120.user_msa.fasta"
    threads:
        config['threads']
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk_align.txt",
        "taxonomy/gtdbtk.log"
    params:
        outdir= "taxonomy"
    shell:
        "GTDBTK_DATA_PATH={GTDBTK_DATA_PATH} ;  "
        "gtdbtk align --identify_dir {params.outdir} --out_dir {params.outdir} "
        "--cpus {threads} &> {log[0]}"


rule classify:
    input:
        rules.align.output,
        genome_dir=genome_dir,
    output:
        "taxonomy/gtdbtk.bac120.summary.tsv",
        scratch_dir=temp(directory(os.path.join(TMPDIR,"gtdb")))
    threads:
        8
    resources:
        mem=250
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk_classify.txt",
        "taxonomy/gtdbtk.log"
    params:
        outdir= "taxonomy",
        extension="fasta",
    shell:
        "GTDBTK_DATA_PATH={GTDBTK_DATA_PATH} ; mkdir {output.scratch_dir} ;"
        "gtdbtk classify --genome_dir {input.genome_dir} --align_dir {params.outdir} "
        "--out_dir {params.outdir} "
        "--extension {params.extension} "
        "--cpus {threads} &> {log[0]}"
#        "--scratch_dir {output.scratch_dir} "

rule infer:
    input:
        "taxonomy/gtdbtk.bac120.user_msa.fasta"
    output:
        "taxonomy/gtdbtk.unrooted.tree"
    threads:
        config['threads']
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk_infer.txt",
        "taxonomy/gtdbtk.log"
    params:
        outdir= "taxonomy"
    shell:
        "GTDBTK_DATA_PATH={GTDBTK_DATA_PATH} ;  "
        "gtdbtk infer --msa_file {input} --out_dir {params.outdir} "
        "--cpus {threads} &> {log[0]}"


# rule download_gtdb:
#     shell:
#         "old_path=$GTDBTK_DATA_PATH ;"
#         " GTDBTK_DATA_PATH={params.db_path} ;"
#         " download-db.sh ;"
#         " ln -s $GTDBTK_DATA_PATH $old_path "
