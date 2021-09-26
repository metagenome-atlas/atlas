### VAMB


localrules:
    combine_contigs,
    vamb,


rule filter_contigs:
    input:
        "{sample}/{sample}_contigs.fasta",
    output:
        temp("Cobinning/filtered_contigs/{sample}.fasta"),
    params:
        min_length=config["cobining_min_contig_length"],
    log:
        "logs/cobinning/filter_contigs/{sample}.log",
    conda:
        "../envs/required_packages.yaml"
    threads: 1
    resources:
        mem=config["simplejob_mem"],
        java_mem=int(int(config["simplejob_mem"] * JAVA_MEM_FRACTION)),
    shell:
        " reformat.sh in={input} "
        " fastaminlen={params.min_length} "
        " out={output} "
        " overwrite=true "
        " -Xmx{resources.java_mem}G 2> {log} "


localrules:
    combine_contigs,


rule combine_contigs:
    input:
        #Trigers rerun if contigs change
        flag=expand("{sample}/{sample}_contigs.fasta", sample=SAMPLES),
        fasta=ancient(expand(rules.filter_contigs.output[0], sample=SAMPLES)),
    output:
        "Cobinning/combined_contigs.fasta.gz",
    log:
        "logs/cobinning/combine_contigs.log",
    params:
        seperator=config["cobinning_separator"],
        samples=SAMPLES,
    threads: 1
    run:
        import gzip as gz

        with gz.open(output[0], "wt") as fout:

            for sample, input_fasta in zip(params.samples, input.fasta):
                with open(input_fasta) as fin:

                    for line in fin:
                        # if line is a header add sample name
                        if line[0] == ">":
                            line = f">{sample}{params.seperator}" + line[1:]
                        # write each line to the combined file
                        fout.write(line)


rule minimap_index:
    input:
        contigs=rules.combine_contigs.output,
    output:
        mmi=temp("Cobinning/combined_contigs.mmi"),
    params:
        index_size="12G",
    resources:
        mem=config["mem"],  # limited num of fatnodes (>200g)
    threads: 3
    log:
        "logs/cobinning/vamb/index.log",
    benchmark:
        "logs/benchmarks/cobinning/mminimap_index.tsv"
    conda:
        "../envs/minimap.yaml"
    shell:
        "minimap2 -I {params.index_size} -t {threads} -d {output} {input} 2> {log}"


rule samtools_dict:
    input:
        contigs=rules.combine_contigs.output,
    output:
        dict="Cobinning/combined_contigs.dict",
    resources:
        mem=config["simplejob_mem"],
        ttime=config["runtime"]["simplejob"],
    threads: 1
    log:
        "logs/cobinning/samtools_dict.log",
    conda:
        "../envs/minimap.yaml"
    shell:
        "samtools dict {input} | cut -f1-3 > {output} 2> {log}"


rule minimap:
    input:
        fq=lambda wildcards: input_paired_only(
            get_quality_controlled_reads(wildcards)
        ),
        mmi="Cobinning/combined_contigs.mmi",
        dict="Cobinning/combined_contigs.dict",
    output:
        bam=temp("Cobinning/mapping/{sample}.unsorted.bam"),
    threads: config["threads"]
    resources:
        mem=config["mem"],
    log:
        "logs/cobinning/mapping/{sample}.minimap.log",
    benchmark:
        "logs/benchmarks/cobinning/mminimap/{sample}.tsv"
    conda:
        "../envs/minimap.yaml"
    shell:
        """minimap2 -t {threads} -ax sr {input.mmi} {input.fq} | grep -v "^@" | cat {input.dict} - | samtools view -F 3584 -b - > {output.bam} 2>{log}"""


ruleorder: sort_bam > minimap > convert_sam_to_bam


rule sort_bam:
    input:
        "Cobinning/mapping/{sample}.unsorted.bam",
    output:
        "Cobinning/mapping/{sample}.sorted.bam",
    params:
        prefix="Cobinning/mapping/tmp.{sample}",
    threads: 2
    resources:
        mem=config["simplejob_mem"],
        time=int(config["runtime"]["simplejob"]),
    log:
        "logs/cobinning/mapping/{sample}.sortbam.log",
    conda:
        "../envs/minimap.yaml"
    shell:
        "samtools sort {input} -T {params.prefix} --threads {threads} -m 3G -o {output} 2>{log}"


rule summarize_bam_contig_depths:
    input:
        bam=expand(rules.sort_bam.output, sample=SAMPLES),
    output:
        "Cobinning/vamb/coverage.jgi.tsv",
    log:
        "logs/cobinning/vamb/combine_coverage.log",
    conda:
        "../envs/metabat.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
    shell:
        "jgi_summarize_bam_contig_depths "
        " --outputDepth {output} "
        " {input.bam} &> {log} "


localrules:
    convert_jgi2vamb_coverage,


rule convert_jgi2vamb_coverage:
    input:
        "Cobinning/vamb/coverage.jgi.tsv",
    output:
        "Cobinning/vamb/coverage.tsv",
    log:
        "logs/cobinning/vamb/convert_jgi2vamb_coverage.log",
    threads: 1
    script:
        "../scripts/convert_jgi2vamb_coverage.py"


rule run_vamb:
    input:
        coverage="Cobinning/vamb/coverage.tsv",
        fasta=rules.combine_contigs.output,
    output:
        directory("Cobinning/vamb/clustering"),
    conda:
        "../envs/vamb.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["long"],
    log:
        "logs/cobinning/vamb/run_vamb.log",
    benchmark:
        "logs/benchmarks/vamb/run_vamb.tsv"
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


localrules: parse_vamb_output

rule parse_vamb_output:
    input:
        rules.run_vamb.output,
    output:
        renamed_clusters="Cobinning/vamb/clusters.tsv.gz",
        cluster_atributions=expand(vamb_cluster_attribution_path, sample=SAMPLES),
    log:
        "logs/cobinning/vamb_parse_output.log",
    params:
        separator=config["cobinning_separator"],
        fasta_extension=".fna",
        output_path=lambda wc: vamb_cluster_attribution_path,
        samples=SAMPLES,
    script:
        "../scripts/parse_vamb.py"


rule vamb:
    input:
        "Cobinning/vamb/clustering",
        "Cobinning/vamb/clusters.tsv.gz",
        expand("{sample}/binning/vamb/cluster_attribution.tsv", sample=SAMPLES),


include: "semibin.smk"
