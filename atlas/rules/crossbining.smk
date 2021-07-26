### VAMB


localrules:
    combine_contigs,
    vamb,


rule vamb:
    input:
        "Crossbinning/vamb/clustering",


rule combine_contigs:
    input:
        expand("{sample}/{sample}_contigs.fasta", sample=SAMPLES),
    output:
        "Crossbinning/combined_contigs.fasta.gz",
    threads: 1
    conda:
        "../envs/vamb.yaml"
    shell:
        "concatenate.py {output} {input} -m 2000 --keepnames"



rule minimap_index:
    input:
        contigs=rules.combine_contigs.output,
    output:
        mmi=temp("Crossbinning/combined_contigs.mmi"),
    params:
        index_size="12G",
    resources:
        mem=config["large_mem"],
    threads: 1
    log:
        "log/vamb/index.log",
    benchmark:
        "log/benchmarks/crossbining/mminimap_index.tsv"
    conda:
        "../envs/minimap.yaml"
    shell:
        "minimap2 -I {params.index_size} -d {output} {input} 2> {log}"


rule samtools_dict:
    input:
        contigs=rules.combine_contigs.output,
    output:
        dict="Crossbinning/mapping/combined_contigs.dict",
    resources:
        mem=config["simplejob_mem"],
    threads: 1
    log:
        "log/crossbining/samtools_dict.log",
    conda:
        "../envs/minimap.yaml"
    shell:
        "samtools dict {input} | cut -f1-3 > {output} 2> {log}"


rule minimap:
    input:
        fq=lambda wildcards: input_paired_only(
            get_quality_controlled_reads(wildcards)
        ),
        mmi="Crossbinning/combined_contigs.mmi",
        dict="Crossbinning/combined_contigs.dict",
    output:
        bam=temp("Crossbinning/mapping/{sample}.bam"),
    threads: config["threads"]
    resources:
        mem = config["mem"],
    log:
        "log/crossbining/mapping/{sample}.minimap.log",
    benchmark:
        "log/benchmarks/crossbining/mminimap/{sample}.tsv"
    conda:
        "../envs/minimap.yaml"
    shell:
        """minimap2 -t {threads} -ax sr {input.mmi} {input.fq} | grep -v "^@" | cat {input.dict} - | samtools view -F 3584 -b - > {output.bam} 2>{log}"""


ruleorder: sort_bam > minimap > convert_sam_to_bam


rule sort_bam:
    input:
        "Crossbinning/mapped/{sample}.bam",
    output:
        temp("Crossbinning/mapped/{sample}.sort.bam"),
    params:
        prefix="Crossbinning/mapped/tmp.{sample}",
    threads: 2
    resources:
        mem=config["simplejob_mem"],
        time=int(config["runtime"]["simple_job"]),
    log:
        "log/crossbining/mapping/{sample}.sortbam.log",
    conda:
        "../envs/minimap.yaml"
    shell:
        "samtools sort {input} -T {params.prefix} --threads {threads} -m 3G -o {output} 2>{log}"


rule summarize_bam_contig_depths:
    input:
        bam=expand(rules.sort_bam.output, sample=SAMPLES),
    output:
        "Crossbinning/vamb/coverage.jgi.tsv",
    log:
        "log/vamb/combine_coverage.log",
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
        "Crossbinning/vamb/coverage.jgi.tsv",
    output:
        "Crossbinning/vamb/coverage.tsv",
    log:
        "log/vamb/convert_jgi2vamb_coverage.log",
    threads: 1
    script:
        "../scripts/convert_jgi2vamb_coverage.py"


rule run_vamb:
    input:
        coverage="Crossbinning/vamb/coverage.tsv",
        fasta=rules.combine_contigs.output,
    output:
        directory("Crossbinning/vamb/clustering"),
    conda:
        "../envs/vamb.yaml"
    threads: 1 #config["threads"]
    resources:
        mem=config["mem"],
        time = config["runtime"]["default"]
    log:
        "log/vamb/run_vamb.log",
    benchmark:
        "log/benchmarks/vamb/run_vamb.tsv"
    params:
        params="-m 2000 --minfasta 500000",
    shell:
        "vamb --outdir {output} "
        " -o '_' "
        " --jgi {input.coverage} "
        " --fasta {input.fasta} "
        " {params.params} "
        "2> {log}"


include: "semibin.smk"
