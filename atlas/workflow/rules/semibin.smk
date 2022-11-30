rule run_semibin:
    input:
        "{sample}/{sample}_contigs.fasta",
        fasta=rules.filter_contigs.output,
        bams=expand(rules.sort_bam.output, sample=SAMPLES),
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
        "SemiBin single_easy_bin "
        " --input-fasta {input.fasta} "
        " --input-bam {input.bams} "
        " --output {params.output_dir} "
        " --threads {threads} "
	" --training-type=self "
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
