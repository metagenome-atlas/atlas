SEMIBIN_DATA_PATH = os.path.join(DBDIR, "SemiBin_GTDB")


localrules:
    semibin_download_gtdb,


rule semibin_download_gtdb:
    output:
        directory(SEMIBIN_DATA_PATH),
    log:
        "logs/download/Semibin.txt",
    conda:
        "../envs/semibin.yaml"
    threads: 1
    shell:
        "SemiBin download_GTDB --reference-db {output} 2> {log}"
        # Semibin 0.2 has the following error https:/github.com/BigDataBiology/SemiBin/issues/31


rule semibin_predict_taxonomy:
    input:
        "{sample}/{sample}_contigs.fasta",
        fasta=rules.filter_contigs.output,
        db=SEMIBIN_DATA_PATH,
    output:
        "Cobinning/SemiBin/{sample}/cannot/cannot_bin.txt",
        "Cobinning/SemiBin/{sample}/mmseqs_contig_annotation/taxonomyResult.tsv",
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem=config["large_mem"],
        time=config["runtime"]["default"],
    log:
        "logs/semibin/predict_taxonomy/{sample}.log",
    benchmark:
        "logs/benchmarks/semibin/predict_taxonomy/{sample}.tsv"
    params:
        output_dir="Cobinning/SemiBin/{sample}",
        name=lambda wc, output: os.path.basename(output[0]).replace(".txt", ""),
        tmp_dir= config['tmpdir']
    shadow:
        "minimal"
    shell:
        " export TMPDIR={params.tmp_dir} &> {log} ;"
        "SemiBin predict_taxonomy "
        " --input-fasta {input.fasta} "
        " --threads {threads} "
        " --output {params.output_dir} "
        " --cannot-name {params.name} "
        " --reference-db {input.db}/GTDB "
        " &>> {log} "


rule semibin_generate_data_multi:
    input:
        fasta=rules.combine_contigs.output,
        bams=expand(rules.sort_bam.output, sample=SAMPLES),
    output:
        expand(
            "Cobinning/SemiBin/samples/{sample}/{files}",
            sample=SAMPLES,
            files=["data.csv", "data_split.csv"],
        ),
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["default"],
    log:
        "logs/semibin/generate_data_multi.log",
    benchmark:
        "logs/benchmarks/semibin/generate_data_multi.tsv"
    params:
        output_dir="Cobinning/SemiBin",
        separator=config["cobinning_separator"],
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
        "{sample}/{sample}_contigs.fasta",
        fasta=rules.filter_contigs.output,
        bams=expand(rules.sort_bam.output, sample=SAMPLES),
        data="Cobinning/SemiBin/samples/{sample}/data.csv",
        data_split="Cobinning/SemiBin/samples/{sample}/data_split.csv",
        cannot_link=rules.semibin_predict_taxonomy.output[0],
    output:
        "Cobinning/SemiBin/{sample}/model.h5",
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["default"],
    log:
        "logs/semibin/train/{sample}.log",
    benchmark:
        "logs/benchmarks/semibin/train/{sample}.tsv"
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
        "{sample}/{sample}_contigs.fasta",
        fasta=rules.filter_contigs.output,
        bams=expand(rules.sort_bam.output, sample=SAMPLES),
        data="Cobinning/SemiBin/samples/{sample}/data.csv",
        model=rules.semibin_train.output[0],
    output:
        directory("Cobinning/SemiBin/{sample}/output_recluster_bins/"),
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["default"],
    log:
        "logs/semibin/bin/{sample}.log",
    benchmark:
        "logs/benchmarks/semibin/bin/{sample}.tsv"
    params:
        output_dir="Cobinning/SemiBin/{sample}/",
        min_bin_kbs=int(config["cobining_min_bin_size"] / 1000),
        extra=config["semibin_options"],
    shell:
        "SemiBin bin "
        " --input-fasta {input.fasta} "
        " --input-bam {input.bams} "
        " --output {params.output_dir} "
        " --threads {threads} "
        " --data {input.data} "
        " --model {input.model} "
        " --minfasta-kbs {params.min_bin_kbs}"
        " {params.extra} "
        " 2> {log}"

localrules: parse_semibin_output

rule parse_semibin_output:
    input:
        rules.run_semibin.output[0],
    output:
        "{sample}/binning/SemiBin/cluster_attribution.tsv",
    log:
        "logs/semibin/parse_output/{sample}.log",
    params:
        extension=".fa",
    script:
        "../scripts/parse_semibin.py"


rule semibin:
    input:
        expand("Cobinning/SemiBin/{sample}/output_recluster_bins/", sample=SAMPLES),
        expand("{sample}/binning/SemiBin/cluster_attribution.tsv", sample=SAMPLES),


# alternative to pretrained model --environment: Environment for the built-in model(human_gut/dog_gut/ocean).”
