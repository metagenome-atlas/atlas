from glob import glob


include: "bin_quality.smk"


rule pileup_for_binning:
    input:
        fasta=get_assembly,
        bam=get_bam,
    output:
        covstats="Intermediate/binning/coverage/{sample}_to_{sample_reads}_coverage_stats.txt",
    params:
        pileup_secondary=(
            "t"
            if config.get("count_multi_mapped_reads", CONTIG_COUNT_MULTI_MAPPED_READS)
            else "f"
        ),
    log:
        "logs/binning/calculate_coverage/pileup_reads_from_{sample_reads}_to_{sample}.log", 
    conda:
        "../envs/required_packages.yaml"
    threads: config["threads"]
    resources:
        mem_mb=config["mem"] * 1000,
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
    shell:
        "pileup.sh "
        " ref={input.fasta} "
        " in={input.bam} "
        " threads={threads} "
        " -Xmx{resources.java_mem}G "
        " covstats={output.covstats} "
        " secondary={params.pileup_secondary} "
        " 2> {log} "


localrules:
    get_contig_coverage_from_bb,
    combine_coverages,


rule get_contig_coverage_from_bb:
    input:
        coverage=rules.pileup_for_binning.output.covstats,
    output:
        temp("Intermediate/binning/contig_coverage/{sample}_to_{sample_reads}.txt"),
    run:
        with open(input[0]) as fi, open(output[0], "w") as fo:
            # header
            next(fi)
            for line in fi:
                toks = line.strip().split("\t")
                print(toks[0], toks[1], sep="\t", file=fo)


rule combine_coverages:
    input:
        covstats=lambda wc: expand(
            "Intermediate/binning/contig_coverage/{sample}_to_{sample_reads}.txt",
            sample_reads=get_alls_samples_of_group(wc),
            sample=wc.sample,
        ),
    output:
        "Intermediate/binning/combined_contig_coverage/{sample}.tsv",
    run:
        from utils.parsers_bbmap import combine_coverages

        combined_cov, _ = combine_coverages(
            input.covstats, get_alls_samples_of_group(wildcards), "Avg_fold"
        )

        combined_cov.T.to_csv(output[0], sep="\t")


## METABAT
rule get_metabat_depth_file:
    input:
        bams=lambda wc: expand(
            get_bam,
            sample_reads=get_alls_samples_of_group(wc),
            sample=wc.sample,
        ),
    output:
        temp("Intermediate/binning/metabat/{sample}/metabat_depth.txt"),
    log:
        "logs/binning/{sample}/binning/jgi_summarize_bam_contig_depths.log",
    conda:
        "../envs/metabat.yaml"
    threads: config["threads"]  # multithreaded trough OMP_NUM_THREADS
    resources:
        mem_mb=config["mem"] * 1000,
    params:
        minid=lambda wc, input: (
            config["cobinning_readmapping_id"] * 100 if len(input.bams) > 1 else 97
        ),
    shell:
        "jgi_summarize_bam_contig_depths "
        " --percentIdentity {params.minid} "
        " --outputDepth {output} "
        " {input.bams} &> {log} "


def get_metabat_sensitivity():
    if config["metabat"]["sensitivity"] == "sensitive":
        return 500
    else:
        200


rule metabat:
    input:
        depth_file=rules.get_metabat_depth_file.output,
        contigs=get_assembly,
    output:
        "Intermediate/binning/metabat/{sample}/cluster_attribution.tmp",
    params:
        sensitivity=get_metabat_sensitivity(),
        min_contig_len=config["metabat"]["min_contig_length"],
        output_prefix="{sample}/binning/bins/bin",
    benchmark:
        "logs/benchmarks/binning/metabat/{sample}.txt"
    log:
        "logs/binning/{sample}/metabat.log",
    conda:
        "%s/metabat.yaml" % CONDAENV
    threads: config["threads"]
    resources:
        mem_mb=config["mem"] * 1000,
    shell:
        """
        metabat2 -i {input.contigs} \
            --abdFile {input.depth_file} \
            --minContig {params.min_contig_len} \
            --numThreads {threads} \
            --maxEdges {params.sensitivity} \
            --saveCls --noBinOut \
            -o {output} \
            &> {log}
        """


rule maxbin:
    input:
        fasta=get_assembly,
        abund="Intermediate/binning/contig_coverage/{sample}_to_{sample}.txt",
    output:
        directory("Intermediate/binning/maxbin/{sample}/intermediate_files"),
        "Intermediate/binning/maxbin/{sample}/{sample}.summary",
        "Intermediate/binning/maxbin/{sample}/{sample}.marker",
        "Intermediate/binning/maxbin/{sample}/{sample}.marker_of_each_bin.tar.gz",
        "Intermediate/binning/maxbin/{sample}/maxbin.log",
    params:
        mi=config["maxbin"]["max_iteration"],
        mcl=config["maxbin"]["min_contig_length"],
        pt=config["maxbin"]["prob_threshold"],
        output_prefix=lambda wc, output: os.path.join(output[0], wc.sample),
    log:
        "logs/binning/{sample}/maxbin.log",
    conda:
        "../envs/maxbin.yaml"
    threads: config["threads"]
    shell:
        """
        mkdir {output[0]} 2> {log}
        run_MaxBin.pl -contig {input.fasta} \
            -abund {input.abund} \
            -out {params.output_prefix} \
            -min_contig_length {params.mcl} \
            -thread {threads} \
            -prob_threshold {params.pt} \
            -max_iteration {params.mi} >> {log}

        mv {params.output_prefix}.summary {output[0]}/.. 2>> {log}
        mv {params.output_prefix}.marker {output[0]}/..  2>> {log}
        mv {params.output_prefix}.marker_of_each_bin.tar.gz {output[0]}/..  2>> {log}
        mv {params.output_prefix}.log {output[0]}/..  2>> {log}

        """


localrules:
    get_bins,


localrules:
    get_unique_cluster_attribution,
    get_maxbin_cluster_attribution,


rule get_unique_cluster_attribution:
    input:
        "Intermediate/binning/{binner}/{sample}/cluster_attribution.tmp",
    output:
        "Binning/cluster_attribution/{binner}/{sample}.tsv",
    run:
        import pandas as pd
        import numpy as np


        d = pd.read_csv(input[0], index_col=0, header=None, sep="\t").squeeze()

        assert (
            type(d) == pd.Series
        ), "expect the input to be a two column file: {}".format(input[0])

        old_cluster_ids = list(d.unique())
        if 0 in old_cluster_ids:
            old_cluster_ids.remove(0)

        map_cluster_ids = dict(
            zip(
                old_cluster_ids,
                utils.gen_names_for_range(
                    len(old_cluster_ids),
                    prefix="{sample}_{binner}_".format(**wildcards),
                ),
            )
        )

        new_d = d.map(map_cluster_ids)
        new_d.dropna(inplace=True)
        if new_d.shape[0] == 0:
            logger.warning(
                f"No bins detected with binner {wildcards.binner} in sample {wildcards.sample}.\n"
                "I add longest contig to make the pipeline continue"
            )

            new_d[f"{wildcards.sample}_0"] = "{sample}_{binner}_1".format(**wildcards)

        new_d.to_csv(output[0], sep="\t", header=False)


#


rule get_maxbin_cluster_attribution:
    input:
        "Intermediate/binning/maxbin/{sample}/intermediate_files",
    output:
        "Intermediate/binning/maxbin/{sample}/cluster_attribution.tmp",
    params:
        file_name=lambda wc, input: "{folder}/{sample}.{{binid}}.fasta".format(
            folder=input[0], **wc
        ),
    run:
        (bin_ids,) = glob_wildcards(params.file_name)
        print("found {} bins".format(len(bin_ids)))
        with open(output[0], "w") as out_file:
            for binid in bin_ids:
                with open(params.file_name.format(binid=binid)) as bin_file:
                    for line in bin_file:
                        if line.startswith(">"):
                            fasta_header = line[1:].strip().split()[0]
                            out_file.write(f"{fasta_header}\t{binid}\n")
                os.remove(params.file_name.format(binid=binid))


rule get_bins:
    input:
        cluster_attribution="Binning/cluster_attribution/{binner}/{sample}.tsv",
        contigs=get_assembly,
    output:
        directory("Intermediate/binning/{binner}/bins/{sample}"),
    conda:
        "../envs/sequence_utils.yaml"
    log:
        "logs/binning/{sample}/get_bins_{binner}.log",
    script:
        "../scripts/get_fasta_of_bins.py"


localrules:
    get_unique_bin_ids,


rule get_unique_bin_ids:
    input:
        "Binning/cluster_attribution/{binner}/{sample}.tsv",
    output:
        "Intermediate/binning/DASTool/{binner}/{sample}.scaffolds2bin",
    shell:
        "cp {input} {output}"


rule run_das_tool:
    input:
        cluster_attribution=expand(
            "Intermediate/binning/DASTool/{binner}/{{sample}}.scaffolds2bin",
            binner=config["binner"],
        ),
        contigs=get_assembly,
        proteins="{sample}/annotation/predicted_genes/{sample}.faa",
    output:
        "Binning/DASTool/{sample}_DASTool_summary.tsv",
        "Binning/DASTool/{sample}_allBins.eval",
        "Binning/DASTool/{sample}_DASTool_contig2bin.tsv",
        cluster_attribution="Binning/cluster_attribution/DASTool/{sample}.tsv",
    threads: config["threads"]
    log:
        "{sample}/logs/binning/DASTool.log",
    conda:
        "%s/DASTool.yaml" % CONDAENV
    params:
        binner_names=",".join(config["binner"]),
        scaffolds2bin=lambda wc, input: ",".join(input.cluster_attribution),
        output_prefix="Binning/DASTool/{sample}",
        score_threshold=config["DASTool"]["score_threshold"],
        megabin_penalty=config["DASTool"]["megabin_penalty"],
        duplicate_penalty=config["DASTool"]["duplicate_penalty"],
    shell:
        " DAS_Tool --outputbasename {params.output_prefix} "
        " --bins {params.scaffolds2bin} "
        " --labels {params.binner_names} "
        " --contigs {input.contigs} "
        " --search_engine diamond "
        " --proteins {input.proteins} "
        " --write_bin_evals "
        " --megabin_penalty {params.megabin_penalty}"
        " --duplicate_penalty {params.duplicate_penalty} "
        " --threads {threads} "
        " --debug "
        " --score_threshold {params.score_threshold} &> {log} "
        " ; cp {params.output_prefix}_DASTool_contig2bin.tsv {output.cluster_attribution} &>> {log}"


#
