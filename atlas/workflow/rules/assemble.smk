





assembly_params = {}

if  (config["data_type"]=="metagenome") and (config.get("assembler", "megahit") == "megahit"):
    assembly_params["megahit"] = {
        "default": "",
        "meta-sensitive": "--presets meta-sensitive",
        "meta-large": " --presets meta-large",
    }
    ASSEMBLY_FRACTIONS = MULTIFILE_FRACTIONS
    if PAIRED_END and config.get("merge_pairs_before_assembly", True):
        ASSEMBLY_FRACTIONS = ["R1", "R2", "me"]

    localrules:
        merge_se_me_for_megahit,

    rule merge_se_me_for_megahit:
        input:
            expand(
                "{{sample}}/assembly/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
                fraction=["se", "me"],
                assembly_preprocessing_steps=assembly_preprocessing_steps,
            ),
        output:
            temp(
                expand(
                    "{{sample}}/assembly/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
                    fraction=["co"],
                    assembly_preprocessing_steps=assembly_preprocessing_steps,
                )
            ),
        shell:
            "cat {input} > {output}"


    def megahit_input_parsing(input):
        Nfiles = len(input)

        if Nfiles == 1:
            out = f"--read {input[0]}"
        else:
            out = f"-1 {input[0]} -2 {input[1]} "

            if Nfiles == 3:
                out += f"--read {input[2]}"
        return out

    rule run_megahit:
        input:
            expand(
                "{{sample}}/assembly/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
                fraction=ASSEMBLY_FRACTIONS,
                assembly_preprocessing_steps=assembly_preprocessing_steps,
            ),
        output:
            temp("{sample}/assembly/megahit/{sample}_prefilter.contigs.fa"),
        benchmark:
            "logs/benchmarks/assembly/megahit/{sample}.txt"
        log:
            "{sample}/logs/assembly/megahit.log",
        params:
            min_count=config.get("megahit_min_count", MEGAHIT_MIN_COUNT),
            k_min=config.get("megahit_k_min", MEGAHIT_K_MIN),
            k_max=config.get("megahit_k_max", MEGAHIT_K_MAX),
            k_step=config.get("megahit_k_step", MEGAHIT_K_STEP),
            merge_level=config.get("megahit_merge_level", MEGAHIT_MERGE_LEVEL),
            prune_level=config.get("megahit_prune_level", MEGAHIT_PRUNE_LEVEL),
            low_local_ratio=config["megahit_low_local_ratio"],
            min_contig_len=config["minimum_contig_length"],
            outdir=lambda wc, output: os.path.dirname(output[0]),
            inputs=lambda wc, input: megahit_input_parsing(input),
            preset=assembly_params["megahit"][config["megahit_preset"]],
        conda:
            "../envs/megahit.yaml"
        threads: config["assembly_threads"]
        resources:
            mem=config["assembly_memory"],
            time=config["runtime"]["assembly"],
        shell:
            """
            rm -r {params.outdir} 2> {log}

            megahit \
            {params.inputs} \
            --tmp-dir {TMPDIR} \
            --num-cpu-threads {threads} \
            --k-min {params.k_min} \
            --k-max {params.k_max} \
            --k-step {params.k_step} \
            --out-dir {params.outdir} \
            --out-prefix {wildcards.sample}_prefilter \
            --min-contig-len {params.min_contig_len} \
            --min-count {params.min_count} \
            --merge-level {params.merge_level} \
            --prune-level {params.prune_level} \
            --low-local-ratio {params.low_local_ratio} \
            --memory {resources.mem}000000000  \
            {params.preset} >> {log} 2>&1
            """

    localrules:
        rename_megahit_output,

    rule rename_megahit_output:
        input:
            "{sample}/assembly/megahit/{sample}_prefilter.contigs.fa",
        output:
            temp("{sample}/assembly/{sample}_raw_contigs.fasta"),
        shell:
            "cp {input} {output}"


else: #spades

    if PAIRED_END:

        ASSEMBLY_FRACTIONS = ["R1", "R2"]
        if config.get("merge_pairs_before_assembly", True):
            ASSEMBLY_FRACTIONS += ["me"]
    else:

        from copy import deepcopy
        ASSEMBLY_FRACTIONS = deepcopy(MULTIFILE_FRACTIONS)

        if config["spades_preset"] == "meta":
            logger.error(
                "Metaspades cannot handle single end libraries. Use another assembler or specify 'spades_preset': normal"
            )
            exit(1)

    assembly_params["spades"] = {"meta": "--meta", "normal": "", "rna": "--rna"}

    def spades_parameters(wc, input):
        if not os.path.exists("{sample}/assembly/params.txt".format(sample=wc.sample)):

            params = {}

            reads = dict(zip(ASSEMBLY_FRACTIONS, input))

            if not PAIRED_END:
                params["inputs"] = " -s {se} ".format(**reads)
            else:
                params["inputs"] = " --pe1-1 {R1} --pe1-2 {R2} ".format(**reads)

                if "se" in ASSEMBLY_FRACTIONS:
                    params["inputs"] += "--pe1-s {se} ".format(**reads)
                if "me" in ASSEMBLY_FRACTIONS:
                    params["inputs"] += "--pe1-m {me} ".format(**reads)

            # Long reads:

            if (config["longread_type"] is not None) & (
                str(config["longread_type"]).lower() != "none"
            ):

                long_read_file = get_files_from_sampleTable(wc.sample, "longreads")[0]
                params["longreads"] = " --{t} {f} ".format(
                    t=config["longread_type"], f=long_read_file
                )
            else:
                params["longreads"] = ""

            params["preset"] = assembly_params["spades"][config["spades_preset"]]
            params["skip_error_correction"] = (
                "--only-assembler" if config["spades_skip_BayesHammer"] else ""
            )
            params["extra"] = config["spades_extra"]

        else:

            params = {
                "inputs": "--restart-from last",
                "preset": "",
                "skip_error_correction": "",
                "extra": "",
                "longreads": "",
            }

        params["outdir"] = "{sample}/assembly".format(sample=wc.sample)

        return params


    rule run_spades:
        input:
            expand(
                "{{sample}}/assembly/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
                fraction=ASSEMBLY_FRACTIONS,
                assembly_preprocessing_steps=assembly_preprocessing_steps,
            ),
        output:
            "{sample}/assembly/contigs.fasta",
            "{sample}/assembly/scaffolds.fasta",
        benchmark:
            "logs/benchmarks/assembly/spades/{sample}.txt"
        params:
            p=lambda wc, input: spades_parameters(wc, input),
            k=config.get("spades_k", SPADES_K),
        log:
            "{sample}/logs/assembly/spades.log",
        conda:
            "../envs/spades.yaml"
        threads: config["assembly_threads"]
        resources:
            mem=config["assembly_memory"],
            time=config["runtime"]["assembly"],
        shell:
            # remove pipeline_state file to create all output files again
            " rm -f {params.p[outdir]}/pipeline_state/stage_*_copy_files 2> {log} ; "
            " "
            "spades.py "
            " --threads {threads} "
            " --memory {resources.mem} "
            " -o {params.p[outdir]} "
            " -k {params.k}"
            " {params.p[preset]} "
            " {params.p[extra]} "
            " {params.p[inputs]} "
            " {params.p[longreads]} "
            " {params.p[skip_error_correction]} "
            " >> {log} 2>&1 "


    use rule run_spades as run_rnaspades with:
        output: 
            "{sample}/assembly/transcripts.fasta"



    localrules:
        rename_spades_output,


    def get_spades_sequences(wc):

        if config["data_type"] == "metatranscriptome":
            sequences= "transcripts"
        else:
            if config["spades_use_scaffolds"]:
                sequences= "scaffolds"
            else:
                sequences= "contigs"

        return "{sample}/assembly/{sequences}.fasta".format(sequences=sequences,**wc)


    rule rename_spades_output:
        input:
            get_spades_sequences
        output:
            temp("{sample}/assembly/{sample}_raw_contigs.fasta"),
        shell:
            "cp {input} {output}"


# standardizes header labels within contig FASTAs
rule rename_contigs:
    input:
        "{sample}/assembly/{sample}_raw_contigs.fasta",
    output:
        "{sample}/assembly/{sample}_prefilter_contigs.fasta",
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads: config.get("simplejob_threads", 1)
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],
    log:
        "{sample}/logs/assembly/post_process/rename_and_filter_size.log",
    params:
        minlength=config["minimum_contig_length"],
    shell:
        "rename.sh "
        " in={input} out={output} ow=t "
        " prefix={wildcards.sample} "
        " minscaf={params.minlength} &> {log} "


rule calculate_contigs_stats:
    input:
        "{sample}/assembly/{sample}_{assembly_step}_contigs.fasta",
    output:
        "{sample}/assembly/contig_stats/{assembly_step}_contig_stats.txt",
    conda:
        "../envs/required_packages.yaml"
    log:
        "{sample}/logs/assembly/post_process/contig_stats_{assembly_step}.log",
    threads: 1
    resources:
        mem=1,
        time=config["runtime"]["simplejob"],
    shell:
        "stats.sh in={input} format=3 out={output} &> {log}"


rule combine_sample_contig_stats:
    input:
        expand(
            "{{sample}}/assembly/contig_stats/{assembly_step}_contig_stats.txt",
            assembly_step=["prefilter", "final"],
        ),
    output:
        "{sample}/assembly/contig_stats.tsv",
    run:
        import os
        import pandas as pd

        c = pd.DataFrame()
        for f in input:
            df = pd.read_csv(f, sep="\t")
            assembly_step = os.path.basename(f).replace("_contig_stats.txt", "")
            c.loc[assembly_step]

        c.to_csv(output[0], sep="\t")


if config["filter_contigs"]:

    ruleorder: align_reads_to_prefilter_contigs > align_reads_to_final_contigs

    rule align_reads_to_prefilter_contigs:
        input:
            query=get_quality_controlled_reads,
            target=rules.rename_contigs.output,
        output:
            bam=temp("{sample}/sequence_alignment/alignment_to_prefilter_contigs.bam"),
        params:
            extra="-x sr",
        log:
            "{sample}/logs/assembly/post_process/align_reads_to_prefiltered_contigs.log",
        threads: config["threads"]
        resources:
            mem_mb=config["mem"] * 1000,
        wrapper:
            "v1.19.0/bio/minimap2/aligner"

    rule pileup_prefilter:
        input:
            fasta="{sample}/assembly/{sample}_prefilter_contigs.fasta",
            bam="{sample}/sequence_alignment/alignment_to_prefilter_contigs.bam",
        output:
            covstats="{sample}/assembly/contig_stats/prefilter_coverage_stats.txt",
        params:
            pileup_secondary="t",
        log:
            "{sample}/logs/assembly/post_process/pilup_prefilter_contigs.log",
        conda:
            "../envs/required_packages.yaml"
        threads: config["threads"]
        resources:
            mem_mb=config["mem"] * 1000,
            java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
        shell:
            "pileup.sh ref={input.fasta} in={input.bam} "
            " threads={threads} "
            " -Xmx{resources.java_mem}G "
            " covstats={output.covstats} "
            " concise=t "
            " secondary={params.pileup_secondary} "
            " 2> {log}"

    rule filter_by_coverage:
        input:
            fasta="{sample}/assembly/{sample}_prefilter_contigs.fasta",
            covstats="{sample}/assembly/contig_stats/prefilter_coverage_stats.txt",
        output:
            fasta="{sample}/assembly/{sample}_final_contigs.fasta",
            removed_names="{sample}/assembly/{sample}_discarded_contigs.fasta",
        params:
            minc=config["minimum_average_coverage"],
            minp=config["minimum_percent_covered_bases"],
            minr=config.get("minimum_mapped_reads", MINIMUM_MAPPED_READS),
            minl=config.get("minimum_contig_length", MINIMUM_CONTIG_LENGTH),
            trim=config.get("contig_trim_bp", CONTIG_TRIM_BP),
        log:
            "{sample}/logs/assembly/post_process/filter_by_coverage.log",
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads: 1
        resources:
            mem=config["simplejob_mem"],
            java_mem=int(config["simplejob_mem"] * JAVA_MEM_FRACTION),
        shell:
            """filterbycoverage.sh in={input.fasta} \
            cov={input.covstats} \
            out={output.fasta} \
            outd={output.removed_names} \
            minc={params.minc} \
            minp={params.minp} \
            minr={params.minr} \
            minl={params.minl} \
            trim={params.trim} \
            -Xmx{resources.java_mem}G 2> {log}"""
# HACK: this makes two copies of the same file



else:  # no filter

    localrules:
        do_not_filter_contigs,

    rule do_not_filter_contigs:
        input:
            rules.rename_contigs.output,
        output:
            "{sample}/assembly/{sample}_final_contigs.fasta",
        threads: 1
        shell:
            "cp {input} {output}"


localrules:
    finalize_contigs,


rule finalize_contigs:
    input:
        "{sample}/assembly/{sample}_final_contigs.fasta",
    output:
        "{sample}/{sample}_contigs.fasta",
    threads: 1
    run:
        os.symlink(os.path.relpath(input[0], os.path.dirname(output[0])), output[0])


# generalized rule so that reads from any "sample" can be aligned to contigs from "sample_contigs"
rule align_reads_to_final_contigs:
    input:
        query=get_quality_controlled_reads,
        target="{sample_contigs}/{sample_contigs}_contigs.fasta",
    output:
        bam="{sample_contigs}/sequence_alignment/{sample}.bam",
    params:
        extra="-x sr",
        sorting="coordinate",
    benchmark:
        "logs/benchmarks/assembly/calculate_coverage/align_reads_to_filtered_contigs/{sample}_to_{sample_contigs}.txt"
    log:
        "{sample_contigs}/logs/assembly/calculate_coverage/align_reads_from_{sample}_to_filtered_contigs.log",
    threads: config["threads"]
    resources:
        mem_mb=config["mem"] * 1000,
    wrapper:
        "v1.19.0/bio/minimap2/aligner"


rule pileup_contigs_sample:
    input:
        fasta="{sample}/{sample}_contigs.fasta",
        bam="{sample}/sequence_alignment/{sample}.bam",
    output:
        covhist="{sample}/assembly/contig_stats/postfilter_coverage_histogram.txt",
        covstats="{sample}/assembly/contig_stats/postfilter_coverage_stats.txt",
        bincov="{sample}/assembly/contig_stats/postfilter_coverage_binned.txt",
    params:
        pileup_secondary=(
            "t"
            if config.get("count_multi_mapped_reads", CONTIG_COUNT_MULTI_MAPPED_READS)
            else "f"
        ),
    benchmark:
        "logs/benchmarks/assembly/calculate_coverage/pileup/{sample}.txt"
    log:
        "{sample}/logs/assembly/calculate_coverage/pilup_final_contigs.log",  # This log file is uesd for report
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads: config.get("threads", 1)
    resources:
        mem=config["mem"],
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
    shell:
        "pileup.sh "
        " ref={input.fasta} "
        " in={input.bam} "
        " threads={threads} "
        " -Xmx{resources.java_mem}G "
        " covstats={output.covstats} "
        " hist={output.covhist} "
        " concise=t "
        " secondary={params.pileup_secondary} "
        " bincov={output.bincov} "
        " 2> {log} "


rule create_bam_index:
    input:
        "{file}.bam",
    output:
        "{file}.bam.bai",
    conda:
        "../envs/required_packages.yaml"
    threads: 1
    resources:
        mem=2 * config["simplejob_threads"],
    shell:
        "samtools index {input}"




rule predict_genes:
    input:
        "{sample}/{sample}_contigs.fasta",
    output:
        fna="{sample}/annotation/predicted_genes/{sample}.fna",
        faa="{sample}/annotation/predicted_genes/{sample}.faa",
        gff="{sample}/annotation/predicted_genes/{sample}.gff",
    conda:
        "%s/prodigal.yaml" % CONDAENV
    log:
        "{sample}/logs/gene_annotation/prodigal.txt",
    benchmark:
        "logs/benchmarks/prodigal/{sample}.txt"
    threads: 1
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],
    shell:
        """
        prodigal -i {input} -o {output.gff} -d {output.fna} \
            -a {output.faa} -p meta -f gff 2> {log}
        """


localrules:
    get_contigs_from_gene_names,


rule get_contigs_from_gene_names:
    input:
        faa="{sample}/annotation/predicted_genes/{sample}.faa",
    output:
        tsv="{sample}/annotation/predicted_genes/{sample}.tsv",
    run:
        header = [
            "gene_id",
            "Contig",
            "Gene_nr",
            "Start",
            "Stop",
            "Strand",
            "Annotation",
        ]
        with open(output.tsv, "w") as tsv:
            tsv.write("\t".join(header) + "\n")
            with open(input.faa) as fin:
                gene_idx = 0
                for line in fin:
                    if line[0] == ">":
                        text = line[1:].strip().split(" # ")
                        old_gene_name = text[0]
                        text.remove(old_gene_name)
                        old_gene_name_split = old_gene_name.split("_")
                        gene_nr = old_gene_name_split[-1]
                        contig_nr = old_gene_name_split[-2]
                        sample = "_".join(
                            old_gene_name_split[: len(old_gene_name_split) - 2]
                        )
                        tsv.write(
                            "{gene_id}\t{sample}_{contig_nr}\t{gene_nr}\t{text}\n".format(
                                text="\t".join(text),
                                gene_id=old_gene_name,
                                i=gene_idx,
                                sample=sample,
                                gene_nr=gene_nr,
                                contig_nr=contig_nr,
                            )
                        )
                        gene_idx += 1
        #



localrules:
    build_assembly_report,
    combine_contig_stats,


rule combine_contig_stats:
    input:
        contig_stats=expand(
            "{sample}/assembly/contig_stats/final_contig_stats.txt", sample=SAMPLES
        ),
        gene_tables=expand(
            "{sample}/annotation/predicted_genes/{sample}.tsv", sample=SAMPLES
        ),
        mapping_logs=expand(
            "{sample}/logs/assembly/calculate_coverage/pilup_final_contigs.log",
            sample=SAMPLES,
        ),
        # mapping logs will be incomplete unless we wait on alignment to finish
        bams=expand("{sample}/sequence_alignment/{sample}.bam", sample=SAMPLES),
    output:
        combined_contig_stats="stats/combined_contig_stats.tsv",
    params:
        samples=SAMPLES,
    log:
        "logs/assembly/combine_contig_stats.log",
    script:
        "../scripts/combine_contig_stats.py"


rule build_assembly_report:
    input:
        combined_contig_stats="stats/combined_contig_stats.tsv",
    output:
        report="reports/assembly_report.html",
    conda:
        "%s/report.yaml" % CONDAENV
    log:
        "logs/assembly/report.log",
    script:
        "../report/assembly_report.py"
