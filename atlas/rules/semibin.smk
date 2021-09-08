SEMIBIN_DATA_PATH = os.path.join(DBDIR, "SemiBin_GTDB")

def get__filtered_sample_contigs(wildcards):

    return "Cobinning/filtered_contigs/{sample}.fasta.gz".format(**wildcards)

localrules: semibin_download_gtdb
rule semibin_download_gtdb:
    output:
        directory(SEMIBIN_DATA_PATH),
    log:
        "log/download/Semibin.txt",
    conda:
        "../envs/semibin.yaml"
    threads: 1
    shell:
        "SemiBin download_GTDB --reference-db {output}/GTDB 2> {log}"

        # Semibin 0.2 has the following error https://github.com/BigDataBiology/SemiBin/issues/31


rule semibin_predict_taxonomy:
    input:
        fasta=get__filtered_sample_contigs,
        db=SEMIBIN_DATA_PATH,
    output:
        "Cobinning/SemiBin/samples/{sample}/cannot/cannot_bin.txt",
        "Cobinning/SemiBin/samples/{sample}/mmseqs_contig_annotation/taxonomyResult.tsv"
    conda:
        "../envs/semibin.yaml"
    threads: 1
    resources:
        mem=config["large_mem"],
        time=config["runtime"]["default"],
    log:
        "log/semibin/predict_taxonomy/{sample}.log",
    benchmark:
        "log/benchmarks/semibin/predict_taxonomy/{sample}.tsv"
    params:
        output_dir= "Cobinning/SemiBin/samples/{sample}",
        name=lambda wc, output: os.path.basename(output[0]).replace('.txt',''),
    shadow:
        "minimal"
    shell:
        "SemiBin predict_taxonomy "
        " --input-fasta {input.fasta} "
        " --output {params.output_dir} "
        " --cannot-name {params.name} "
        " --reference-db {input.db}/GTDB "
        " > {log} 2> {log}"


rule semibin_generate_data_multi:
    input:
        fasta=rules.combine_contigs.output,
        bams=expand(rules.sort_bam.output, sample=SAMPLES),
    output:
        expand("Cobinning/SemiBin/samples/{sample}/{files}",
        sample=SAMPLES,
        files=  ["data.csv","data_split.csv"]
        )
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["default"],
    log:
        "log/semibin/generate_data_multi.log",
    benchmark:
        "log/benchmarks/semibin/generate_data_multi.tsv"
    params:
        output_dir="Cobinning/SemiBin",
        separator=config['cobinning_separator'],
    shell:
        "SemiBin generate_data_multi "
        " --input-fasta {input.fasta} "
        " --input-bam {input.bams} "
        " --output {params.output_dir} "
        " --threads {threads} "
        " --separator {params.separator} "
        " 2> {log}"


rule semibin_train:
    input:
        fasta= get__filtered_sample_contigs,
        bams=expand(rules.sort_bam.output, sample=SAMPLES),
        data="Cobinning/SemiBin/samples/{sample}/data.csv",
        data_split="Cobinning/SemiBin/samples/{sample}/data_split.csv",
        cannot_link=rules.semibin_predict_taxonomy.output[0],
    output:
        "Cobinning/SemiBin/samples/{sample}/model.h5",
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["default"],
    log:
        "log/semibin/train/{sample}.log",
    benchmark:
        "log/benchmarks/semibin/train/{sample}.tsv"
    params:
        output_dir=lambda wc, output: os.path.dirname(output[0]),
        extra=" --epoches 20 --mode single ",
    shell:
        "SemiBin train "
        " --input-fasta {input.fasta} "
        " --input-bam {input.bams} "
        " --output {params.output_dir} "
        " --threads {threads} "
        " --data {input.data} "
        " --data-split {input.data_split} "
        " --cannot-link {input.cannot_link} "
        " --mode single "
        " {params.extra} "
        " 2> {log}"


rule run_semibin:
    input:
        fasta=get__filtered_sample_contigs,
        bams=expand(rules.sort_bam.output, sample=SAMPLES),
        data="Cobinning/SemiBin/samples/{sample}/data.csv",
        model=rules.semibin_train.output[0],
    output:
        touch("Cobinning/SemiBin/samples/{sample}/finished"),
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["default"],
    log:
        "log/semibin/bin/{sample}.log",
    benchmark:
        "log/benchmarks/semibin/bin/{sample}.tsv"
    params:
        output_dir="Cobinning/SemiBin",
        min_bin_kbs= config["cobining_min_bin_size"] // 1000,
        extra= config['semibin_options'],
    shell:
        "SemiBin train "
        " --input-fasta {input.fasta} "
        " --input-bam {input.bams} "
        " --output {params.output_dir} "
        " --threads {threads} "
        " --data {input.data} "
        " --model {input.model} "
        " --minfasta-kbs {params.min_bin_kbs}"
        " {params.extra} "
        " 2> {log}"


rule semibin:
    input:
        expand("Cobinning/SemiBin/samples/{sample}/finished", sample=SAMPLES)
# alternative to pretrained model --environment: Environment for the built-in model(human_gut/dog_gut/ocean).”
