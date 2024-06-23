
rule semibin_generate_data_multi:
    input:
        fasta=rules.combine_contigs.output,
        bams=get_bams_of_bingroup,
    output:
        directory("Intermediate/cobinning/{bingroup}/semibin/data_multi"),
        # expand(
        #     "Cobinning/SemiBin/samples/{sample}/{files}",
        #     sample=SAMPLES,
        #     files=["data.csv", "data_split.csv"],
        # ),
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem_mb=config["mem"] * 1000,
        time_min=60 * config["runtime"]["default"],
    log:
        "logs/semibin/{bingroup}/generate_data_multi.log",
    benchmark:
        "logs/benchmarks/semibin/{bingroup}/generate_data_multi.tsv"
    params:
        # output_dir="Cobinning/SemiBin",
        separator=config["cobinning_separator"],
    shell:
        "SemiBin generate_sequence_features_multi"
        " --input-fasta {input.fasta} "
        " --input-bam {input.bams} "
        " --output {output} "
        " --threads {threads} "
        " --separator {params.separator} "
        " 2> {log}"


rule semibin_train:
    input:
        flag=get_assembly,
        fasta_sample=rules.filter_contigs.output[0],
        bams=get_bams_of_bingroup,
        data_folder=rules.semibin_generate_data_multi.output[0],
    output:
        "Intermediate/cobinning/{bingroup}/semibin/models/{sample}/model.h5",
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem_mb=config["mem"] * 1000,
        time_min=60 * config["runtime"]["default"],
    log:
        "logs/semibin/{bingroup}/train/{sample}.log",
    benchmark:
        "logs/benchmarks/semibin/{bingroup}/train/{sample}.tsv"
    params:
        output_dir=lambda wc, output: os.path.dirname(output[0]),
        data=lambda wc, input: Path(input.data_folder)
        / "samples"
        / wc.sample
        / "data.csv",
        data_split=lambda wc, input: Path(input.data_folder)
        / "samples"
        / wc.sample
        / "data_split.csv",
        extra=config["semibin_train_extra"],
    shell:
        "SemiBin train_self "
        " --output {params.output_dir} "
        " --threads {threads} "
        " --data {params.data} "
        " --data-split {params.data_split} "
        " {params.extra} "
        " 2> {log}"


def semibin_input(wildcards):
    bingroup_of_sample = sampleTable.loc[wildcards.sample, "BinGroup"]
    samples_of_bingroup = sampleTable.query(
        f'BinGroup=="{bingroup_of_sample}"'
    ).index.tolist()

    assert len(samples_of_bingroup) > 1

    mapping = dict(
        fasta=rules.filter_contigs.output[0].format(**wildcards),
        bams=expand(
            "Intermediate/cobinning/{bingroup}/bams/{sample}.sorted.bam",
            sample=samples_of_bingroup,
            bingroup=bingroup_of_sample,
        ),
        data_folder=rules.semibin_generate_data_multi.output[0].format(
            bingroup=bingroup_of_sample, **wildcards
        ),
        model=rules.semibin_train.output[0].format(
            bingroup=bingroup_of_sample, **wildcards
        ),
    )

    return mapping


rule run_semibin:
    input:
        unpack(semibin_input),
    output:
        # contains no info to bingroup
        directory(
            "Intermediate/cobinning/semibin_output/{sample}/output_recluster_bins/"
        ),
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem_mb=config["mem"] * 1000,
        time_min=60 * config["runtime"]["default"],
    log:
        "logs/semibin/bin/{sample}.log",
    benchmark:
        "logs/benchmarks/semibin/bin/{sample}.tsv"
    params:
        output_dir=lambda wc, output: os.path.dirname(output[0]),
        data=lambda wc, input: Path(input.data_folder)
        / "samples"
        / wc.sample
        / "data.csv",
        min_bin_kbs=int(config["cobining_min_bin_size"] / 1000),
        extra=config["semibin_options"],
    shell:
        "SemiBin bin "
        " --input-fasta {input.fasta} "
        " --output {params.output_dir} "
        " --threads {threads} "
        " --data {params.data} "
        " --model {input.model} "
        " --minfasta-kbs {params.min_bin_kbs}"
        " {params.extra} "
        " 2> {log}"


localrules:
    parse_semibin_output,


ruleorder: parse_semibin_output > get_unique_cluster_attribution


rule parse_semibin_output:
    input:
        rules.run_semibin.output[0],
    output:
        "{sample}/binning/SemiBin/cluster_attribution.tsv",
    conda:
        "../envs/semibin.yaml"
    log:
        "logs/semibin/parse_output/{sample}.log",
    params:
        extension=".fa",
    script:
        "../scripts/parse_semibin.py"


rule semibin:
    input:
        expand("{sample}/binning/SemiBin/cluster_attribution.tsv", sample=SAMPLES),
