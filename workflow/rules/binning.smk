from glob import glob

BINNING_CONTIGS = "{sample}/{sample}_contigs.fasta"


include: "bin_quality.smk"


rule pileup_for_binning:
    input:
        fasta=BINNING_CONTIGS,
        bam="{sample}/sequence_alignment/{sample_reads}.bam",
    output:
        covstats="{sample}/binning/coverage/{sample_reads}_coverage_stats.txt",
    params:
        pileup_secondary=(
            "t"
            if config.get("count_multi_mapped_reads", CONTIG_COUNT_MULTI_MAPPED_READS)
            else "f"
        ),
        minmapq=config["minimum_map_quality"]
    log:
        "{sample}/logs/binning/calculate_coverage/pileup_reads_from_{sample_reads}_to_filtered_contigs.log",  # this file is udes for assembly report
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
        coverage="{sample}/binning/coverage/{sample_reads}_coverage_stats.txt",
    output:
        temp("{sample}/binning/coverage/{sample_reads}_coverage.txt"),
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
            "{sample}/binning/coverage/{sample_reads}_coverage_stats.txt",
            sample_reads=get_alls_samples_of_group(wc),
            sample=wc.sample,
        ),
    output:
        "{sample}/binning/coverage/combined_coverage.tsv",
    run:
        from utils.parsers_bbmap import combine_coverages

        combined_cov, _ = combine_coverages(
            input.covstats, get_alls_samples_of_group(wildcards), "Avg_fold"
        )

        combined_cov.T.to_csv(output[0], sep="\t")


## CONCOCT
rule run_concoct:
    input:
        coverage="{sample}/binning/coverage/combined_coverage.tsv",
        fasta=BINNING_CONTIGS,
    output:
        "{{sample}}/binning/concoct/intermediate_files/clustering_gt{}.csv".format(
            config["concoct"]["min_contig_length"]
        ),
    params:
        basename=lambda wc, output: os.path.dirname(output[0]),
        Nexpected_clusters=config["concoct"]["Nexpected_clusters"],
        read_length=config["concoct"]["read_length"],
        min_length=config["concoct"]["min_contig_length"],
        niterations=config["concoct"]["Niterations"],
    log:
        "{sample}/binning/concoct/intermediate_files/log.txt",
    conda:
        "%s/concoct.yaml" % CONDAENV
    threads: 10  # concoct uses 10 threads by default, wit for update: https://github.com/BinPro/CONCOCT/issues/177
    resources:
        mem=config["mem"],
    shell:
        """
        concoct -c {params.Nexpected_clusters} \
            --coverage_file {input.coverage} \
            --composition_file {input.fasta} \
            --basename {params.basename} \
            --read_length {params.read_length} \
            --length_threshold {params.min_length} \
            --converge_out \
            --iterations {params.niterations}
        """


localrules:
    convert_concoct_csv_to_tsv,


rule convert_concoct_csv_to_tsv:
    input:
        rules.run_concoct.output[0],
    output:
        "{sample}/binning/concoct/cluster_attribution.tmp",
    run:
        with open(input[0]) as fin, open(output[0], "w") as fout:
            for line in fin:
                fout.write(line.replace(",", "\t"))


## METABAT
rule get_metabat_depth_file:
    input:
        bam=lambda wc: expand(
            "{sample}/sequence_alignment/{sample_reads}.bam",
            sample_reads=get_alls_samples_of_group(wc),
            sample=wc.sample,
        ),
    output:
        temp("{sample}/binning/metabat/metabat_depth.txt"),
    log:
        "{sample}/binning/metabat/metabat.log",
    conda:
        "%s/metabat.yaml" % CONDAENV
    threads: config["threads"]
    resources:
        mem=config["mem"],
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output} {input.bam} \
            &> {log}
        """


def get_metabat_sensitivity():
    if config["metabat"]["sensitivity"] == "sensitive":
        return 500
    else:
        200


rule metabat:
    input:
        depth_file=rules.get_metabat_depth_file.output,
        contigs=BINNING_CONTIGS,
    output:
        "{sample}/binning/metabat/cluster_attribution.tmp",
    params:
        sensitivity=get_metabat_sensitivity(),
        min_contig_len=config["metabat"]["min_contig_length"],
        output_prefix="{sample}/binning/bins/bin",
    benchmark:
        "logs/benchmarks/binning/metabat/{sample}.txt"
    log:
        "{sample}/logs/binning/metabat.txt",
    conda:
        "%s/metabat.yaml" % CONDAENV
    threads: config["threads"]
    resources:
        mem=config["mem"],
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
        fasta=BINNING_CONTIGS,
        abund="{sample}/binning/coverage/{sample}_coverage.txt",
    output:
        directory("{sample}/binning/maxbin/intermediate_files"),
        "{sample}/binning/maxbin/{sample}.summary",
        "{sample}/binning/maxbin/{sample}.marker",
        "{sample}/binning/maxbin/{sample}.marker_of_each_bin.tar.gz",
        "{sample}/binning/maxbin/{sample}.log",
    params:
        mi=config["maxbin"]["max_iteration"],
        mcl=config["maxbin"]["min_contig_length"],
        pt=config["maxbin"]["prob_threshold"],
        output_prefix=lambda wc, output: os.path.join(output[0], wc.sample),
    log:
        "{sample}/logs/binning/maxbin.log",
    conda:
        "%s/maxbin.yaml" % CONDAENV
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
        "{sample}/binning/{binner}/cluster_attribution.tmp",
    output:
        "{sample}/binning/{binner}/cluster_attribution.tsv",
    run:
        import pandas as pd
        import numpy as np


        d = pd.read_csv(input[0], index_col=0, squeeze=True, header=None, sep="\t")

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
                "I add longest contig to make the pipline continue"
            )

            new_d[f"{wildcards.sample}_0"] = "{sample}_{binner}_1".format(**wildcards)

        new_d.to_csv(output[0], sep="\t", header=False)


#


rule get_maxbin_cluster_attribution:
    input:
        "{sample}/binning/maxbin/intermediate_files",
    output:
        "{sample}/binning/maxbin/cluster_attribution.tmp",
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
        cluster_attribution="{sample}/binning/{binner}/cluster_attribution.tsv",
        contigs=BINNING_CONTIGS,
    output:
        directory("{sample}/binning/{binner}/bins"),
    conda:
        "../envs/sequence_utils.yaml"
    log:
        "{sample}/logs/binning/get_bins_{binner}.log",
    script:
        "get_fasta_of_bins.py"


localrules:
    get_unique_bin_ids,


rule get_unique_bin_ids:
    input:
        "{sample}/binning/{binner}/cluster_attribution.tsv",
    output:
        "{sample}/binning/DASTool/{binner}.scaffolds2bin",
    shell:
        "cp {input} {output}"


rule run_das_tool:
    input:
        cluster_attribution=expand(
            "{{sample}}/binning/DASTool/{binner}.scaffolds2bin",
            binner=config["binner"],
        ),
        contigs=BINNING_CONTIGS,
        proteins="{sample}/annotation/predicted_genes/{sample}.faa",
    output:
        "{sample}/binning/DASTool/{sample}_DASTool_summary.tsv",
        "{sample}/binning/DASTool/{sample}_allBins.eval",
        cluster_attribution="{sample}/binning/DASTool/cluster_attribution.tsv",
    threads: config["threads"]
    log:
        "{sample}/logs/binning/DASTool.log",
    conda:
        "%s/DASTool.yaml" % CONDAENV
    params:
        binner_names=",".join(config["binner"]),
        scaffolds2bin=lambda wc, input: ",".join(input.cluster_attribution),
        output_prefix="{sample}/binning/DASTool/{sample}",
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
        " ; mv {params.output_prefix}_DASTool_contig2bin.tsv {output.cluster_attribution} &>> {log}"


#
