
rule semibin_generate_data_multi:
    input:
        fasta=rules.combine_contigs.output,
        bams=lambda wc: expand(rules.sort_bam.output, sample= get_samples_of_bingroup(wc)),
    output:
        directory("Intermediate/cobinning/{bingroup}/semibin/data_multi")
        # expand(
        #     "Cobinning/SemiBin/samples/{sample}/{files}",
        #     sample=SAMPLES,
        #     files=["data.csv", "data_split.csv"],
        # ),
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["default"],
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
        flag = "{sample}/{sample}_contigs.fasta",
        fasta_sample = rules.filter_contigs.output,
        bams= rules.semibin_generate_data_multi.input.bams,
        data_folder= rules.semibin_generate_data_multi.output[0],
    output:
        "Intermediate/cobinning/{bingroup}/semibin/models/{sample}/model.h5",
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["default"],
    log:
        "logs/semibin/{bingroup}/train/{sample}.log",
    benchmark:
        "logs/benchmarks/semibin/{bingroup}/train/{sample}.tsv"
    params:
        output_dir=lambda wc, output: os.path.dirname(output[0]),
        data = lambda wc, input: Path(input.data_folder)/"samples"/wc.sample/"data.csv",
        data_split = lambda wc, input: Path(input.data_folder)/"samples"/wc.sample/"data_split.csv",
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

    bingroup_of_sample = sampleTable.loc[wildcards.sample, "bingroup"]
    samples_of_bingroup = sampleTable.query(f'BinGroup=="{bingroup_of_sample}"').index.tolist()


    return dict(
        flag= "{sample}/{sample}_contigs.fasta",
        fasta = rules.filter_contigs.output,
        bams =lambda wc: expand(rules.sort_bam.output, sample= samples_of_bingroup),
        data_folder = rules.semibin_generate_data_multi.output[0].format(bingroup=bingroup_of_sample),
        model = rules.semibin_train.output[0].format(bingroup=bingroup_of_sample, sample=wildcards.sample),
    )

rule run_semibin:
    input:
        unpack(semibin_input),
    output:
        # contains no info to bingroup
        directory("Intermediate/cobinning/semibin_output/{sample}/output_recluster_bins/"),
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
        output_dir= lambda wc, output: os.path.dirname(output[0])
        data = lambda wc, input: Path(input.data_folder)/"samples"/wc.sample/"data.csv",
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
