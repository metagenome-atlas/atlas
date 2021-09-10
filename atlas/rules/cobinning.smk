### VAMB


localrules:
    combine_contigs,
    vamb,


rule vamb:
    input:
        "Cobinning/vamb/clustering",

# def should_add_seperator():
#
#     seperator = config['cobinning_seperator']
#
#     if any('_' in SAMPLES) and (separator=='_'):
#
#         logger.error("You are trying to do cobinning. \n"
#                      "The Samplenames contain '_' which might lead to confusion. \n"
#                      "Add the folowing option to the config.yaml:\n\n"
#                      "   cobinning_separator: ':' "
#                      )
#
#         exit(1)
#
#
#     return  "f" if seperator=='_' else 't'



rule filter_contigs:
    input:
        "{sample}/{sample}_contigs.fasta"
    output:
        "Cobinning/filtered_contigs/{sample}.fasta.gz"

    params:
        min_length= config['cobining_min_contig_length'],
        prefix= lambda wc: wc.Sample + config['cobinning_separator'] ,
        addprefix = "f" if config['cobinning_separator']=='_' else 't'
    log:
        "logs/cobinning/filter_contigs/{sample}.log"
    conda:
        "../envs/required_packages.yaml"
    threads: 1
    resources:
        mem=config["simplejob_mem"],
        java_mem=int(int(config["simplejob_mem"] * JAVA_MEM_FRACTION)),
    shell:
        " rename.sh in={input} "
        " prefix={params.prefix} "
        " addprefix={params.addprefix} "
        " minscaf={params.min_length} "
        " out={output} "
        " overwrite=true "
        " -Xmx{resources.java_mem}G 2> {log} "




localrules: combine_contigs
rule combine_contigs:
    input:
        ancient(expand("Cobinning/filtered_contigs/{sample}.fasta.gz", sample=SAMPLES)),
    output:
        "Cobinning/combined_contigs.fasta.gz",
    log:
        "logs/cobinning/combine_contigs.log",
    threads: 1
    run:
        from utils.io import cat_files
        cat_files(input, output[0])


rule minimap_index:
    input:
        contigs=rules.combine_contigs.output,
    output:
        mmi=temp("Cobinning/combined_contigs.mmi"),
    params:
        index_size="12G",
    resources:
        mem=config["mem"], # limited num of fatnodes (>200g)
    threads: 1
    log:
        "logs/cobinning/vamb/index.log",
    benchmark:
        "logs/benchmarks/cobinning/mminimap_index.tsv"
    conda:
        "../envs/minimap.yaml"
    shell:
        "minimap2 -I {params.index_size} -d {output} {input} 2> {log}"


rule samtools_dict:
    input:
        contigs=rules.combine_contigs.output,
    output:
        dict="Cobinning/combined_contigs.dict",
    resources:
        mem=config["simplejob_mem"],
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
        time=int(config["runtime"]["simple_job"]),
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
        params="-m 2000 ", #TODO: I don't know what this is for
        minfasta=config["cobining_min_contig_length"],
        separator= config['cobinning_separator']
    shell:
        "vamb --outdir {output} "
        " --minfasta {params.minfasta} "
        " -o '{params.separator}' "
        " --jgi {input.coverage} "
        " --fasta {input.fasta} "
        " {params.params} "
        "2> {log}"


include: "semibin.smk"
