### VAMB


localrules:
    combine_contigs,
    vamb,


rule filter_contigs:
    input:
        get_assembly,
    output:
        "Intermediate/cobinning/filtered_contigs/{sample}.fasta.gz",
    params:
        min_length=config["cobining_min_contig_length"],
    log:
        "logs/cobinning/filter_contigs/{sample}.log",
    conda:
        "../envs/required_packages.yaml"
    threads: 1
    resources:
        mem=config["simplejob_mem"],
        java_mem=int(config["simplejob_mem"] * JAVA_MEM_FRACTION),
    shell:
        " reformat.sh in={input} "
        " fastaminlen={params.min_length} "
        " out={output} "
        " overwrite=true "
        " -Xmx{resources.java_mem}G 2> {log} "


def get_samples_of_bingroup(wildcards):

    samples_of_group= sampleTable.query(f'BinGroup=="{wildcards.bingroup}"').index.tolist()

    return samples_of_group


def get_filtered_contigs_of_bingroup(wildcards):

    samples_of_group = get_samples_of_bingroup(wildcards)

    if len(samples_of_group) <= 5:
        raise ValueError(
            f"Bin group {wildcards.bingroup} has {len(samples_of_group)} less than 5 samples."
            "For cobinning we reccomend at least 5 samples per bin group."
            "Adapt the sample.tsv to set BinGroup of size [5- 1000]"
        )

    return expand(rules.filter_contigs.output[0], sample=samples_of_group)


def get_bams_of_bingroup(wildcards):

    samples_of_group = get_samples_of_bingroup(wildcards)

    return expand(
        "Intermediate/cobinning/{bingroup}/bams/{sample}.sorted.bam",
        sample=samples_of_group,
        bingroup=wildcards.bingroup,
    )


rule combine_contigs:
    input:
        fasta= get_filtered_contigs_of_bingroup,
    output:
        "Intermediate/cobinning/{bingroup}/combined_contigs.fasta.gz",
    log:
        "logs/cobinning/{bingroup}/combine_contigs.log",
    params:
        seperator=config["cobinning_separator"],
        samples=get_samples_of_bingroup,
    threads: 1
    run:
        import gzip as gz

        with gz.open(output[0], "wb") as fout:
            for sample, input_fasta in zip(params.samples, input.fasta):
                with gz.open(input_fasta, "rb") as fin:
                    for line in fin:
                        # if line is a header add sample name
                        if line[0] == ord('>'):
                            line = f">{sample}{params.seperator}".encode() + line[1:]
                        # write each line to the combined file
                        fout.write(line)


rule minimap_index:
    input:
        contigs=rules.combine_contigs.output,
    output:
        mmi=temp("Intermediate/cobinning/{bingroup}/combined_contigs.mmi"),
    params:
        index_size="12G",
    resources:
        mem=config["mem"],  # limited num of fatnodes (>200g)
    threads: config["simplejob_threads"],
    log:
        "logs/cobinning/{bingroup}/minimap_index.log",
    benchmark:
        "logs/benchmarks/cobinning/{bingroup}/mminimap_index.tsv"
    conda:
        "../envs/minimap.yaml"
    shell:
        "minimap2 -I {params.index_size} -t {threads} -d {output} {input} 2> {log}"


rule samtools_dict:
    input:
        contigs=rules.combine_contigs.output,
    output:
        dict="Intermediate/cobinning/{bingroup}/combined_contigs.dict",
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],
    threads: 1
    log:
        "logs/cobinning/{bingroup}/samtools_dict.log",
    conda:
        "../envs/minimap.yaml"
    shell:
        "samtools dict {input} | cut -f1-3 > {output} 2> {log}"


rule minimap:
    input:
        fq=get_quality_controlled_reads,
        mmi="Intermediate/cobinning/{bingroup}/combined_contigs.mmi",
        dict="Intermediate/cobinning/{bingroup}/combined_contigs.dict",
    output:
        bam=temp("Intermediate/cobinning/{bingroup}/bams/{sample}.unsorted.bam"),
    threads: config["threads"]
    resources:
        mem=config["mem"],
    log:
        "logs/cobinning/{bingroup}/mapping/minimap/{sample}.log",
    benchmark:
        "logs/benchmarks/cobinning/{bingroup}/mminimap/{sample}.tsv"
    conda:
        "../envs/minimap.yaml"
    shell:
        """minimap2 -t {threads} -ax sr {input.mmi} {input.fq} | grep -v "^@" | cat {input.dict} - | samtools view -F 3584 -b - > {output.bam} 2>{log}"""

    # samtools filters out secondary alignments

ruleorder: sort_bam > minimap


rule sort_bam:
    input:
        "Intermediate/cobinning/{bingroup}/bams/{sample}.unsorted.bam",
    output:
        "Intermediate/cobinning/{bingroup}/bams/{sample}.sorted.bam",
    params:
        prefix="Intermediate/cobinning/{bingroup}/bams/tmp.{sample}",
    threads: 2
    resources:
        mem_mb=config["simplejob_mem"] *1000,
        time=int(config["runtime"]["simplejob"]),
    log:
        "logs/cobinning/{bingroup}/mapping/sortbam/{sample}.log",
    conda:
        "../envs/minimap.yaml"
    shell:
        "samtools sort {input} -T {params.prefix} --threads {threads} -m 3G -o {output} 2>{log}"


rule summarize_bam_contig_depths:
    input:
        bams=get_bams_of_bingroup,
    output:
        "Intermediate/cobinning/{bingroup}/coverage.jgi.tsv",
    log:
        "logs/cobinning/{bingroup}/combine_coverage.log",
    conda:
        "../envs/metabat.yaml"
    threads: 1
    benchmark:
        "logs/benchmarks/cobinning/{bingroup}/summarize_bam_contig_depths.tsv"
    resources:
        mem_mb=config["mem"]*1000,
    shell:
        "jgi_summarize_bam_contig_depths "
        " --outputDepth {output} "
        " {input.bams} &> {log} "


localrules:
    convert_jgi2vamb_coverage,


rule convert_jgi2vamb_coverage:
    input:
        rules.summarize_bam_contig_depths.output,
    output:
        "Intermediate/cobinning/{bingroup}/coverage.tsv",
    log:
        "logs/cobinning/{bingroup}/convert_jgi2vamb_coverage.log",
    benchmark:
        "logs/benchmarks/cobinning/{bingroup}/convert_jgi2vamb_coverage.tsv"
    threads: 1
    script:
        "../scripts/convert_jgi2vamb_coverage.py"


rule run_vamb:
    input:
        coverage="Intermediate/cobinning/{bingroup}/coverage.tsv",
        fasta=rules.combine_contigs.output,
    output:
        directory("Intermediate/cobinning/{bingroup}/vamb_output"),
    conda:
        "../envs/vamb.yaml"
    threads: config["threads"]
    resources:
        mem_mb=config["mem"]*1000,
        time=config["runtime"]["long"],
    log:
        "logs/cobinning/run_vamb/{bingroup}.log",
    benchmark:
        "logs/benchmarks/run_vamb/{bingroup}.tsv"
    params:
        mincontig=config["cobining_min_contig_length"],  # min contig length for binning
        minfasta=config["cobining_min_bin_size"],  # min bin size for output
        separator=config["cobinning_separator"],
    shell:
        "vamb --outdir {output} "
        " -m {params.mincontig} "
        " --minfasta {params.minfasta} "
        " -o '{params.separator}' "
        " --jgi {input.coverage} "
        " --fasta {input.fasta} "
        "2> {log}"


vamb_cluster_attribution_path = "{sample}/binning/vamb/cluster_attribution.tsv"


localrules:
    parse_vamb_output,


rule parse_vamb_output:
    input:
        expand(rules.run_vamb.output, bingroup=sampleTable.BinGroup.unique()),
    output:
        renamed_clusters="Binning/vamb/vamb_clusters.tsv.gz",
        cluster_atributions=expand(vamb_cluster_attribution_path, sample=SAMPLES),
    log:
        "logs/cobinning/vamb_parse_output.log",
    params:
        separator=config["cobinning_separator"],
        fasta_extension=".fna",
        output_path=lambda wc: vamb_cluster_attribution_path,  # path with {sample} to replace
        samples=SAMPLES,
        bingroups = sampleTable.BinGroup.unique()
    conda:
        "../envs/fasta.yaml"
    script:
        "../scripts/parse_vamb.py"


rule vamb:
    input:
        "Binning/vamb/vamb_clusters.tsv.gz",
        expand("{sample}/binning/vamb/cluster_attribution.tsv", sample=SAMPLES),


include: "semibin.smk"
