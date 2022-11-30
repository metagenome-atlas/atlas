
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
        "SemiBin generate_sequence_features_multi"
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
        extra=" --epochs 20 --mode single ",
    shell:
        "SemiBin train_self "
        " --output {params.output_dir} "
        " --threads {threads} "
        " --data {input.data} "
        " --data-split {input.data_split} "
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
        " --output {params.output_dir} "
        " --threads {threads} "
        " --data {input.data} "
        " --model {input.model} "
        " --minfasta-kbs {params.min_bin_kbs}"
        " {params.extra} "
        " 2> {log}"


localrules:
    parse_semibin_output,


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
