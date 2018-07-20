

localrules: get_contig_covarage_from_bb
rule get_contig_covarage_from_bb:
    input:
        coverage = BB_COVERAGE_FILE
    output:
        "{sample}/binning/contig_coverage.tsv",
    run:
        with open(input[0]) as fi, open(output[0], "w") as fo:
            # header
            next(fi)
            for line in fi:
                toks = line.strip().split("\t")
                print(toks[0], toks[1], sep="\t", file=fo)


## CONCOCT


rule run_concoct:
    input:
        coverage = "{sample}/binning/contig_coverage.tsv",
        fasta = BINNING_CONTIGS
    output:
        "{{sample}}/binning/concoct/intermediate_files/clustering_gt{}.csv".format(config["binning_min_contig_length"])
    params:
        basename= lambda wc, output: os.path.dirname(output[0]),
        Nexpected_clusters= config['concoct']['Nexpected_clusters'],
        read_length= config['concoct']['read_length'],
        min_length=config["binning_min_contig_length"],
        niterations=config["concoct"]["Niterations"]
    log:
        "{sample}/logs/binning/concoct/log.txt"
    conda:
        "%s/concoct.yaml" % CONDAENV
    threads:
        10 # concoct uses 10 threads by default, wit for update: https://github.com/BinPro/CONCOCT/issues/177
    resources:
        mem = config["java_mem"]
    shell:
        """
        concoct -c {params.Nexpected_clusters}\
            --coverage_file {input.coverage}\
            --composition_file {input.fasta}\
            --basename {params.basename}\
            --read_length {params.read_length} \
            --length_threshold {params.min_length}\
            --converge_out \
            --iterations {params.niterations}
        """
localrules: get_cluster_attribution_concoct

rule concoct:
    input:
        rules.run_concoct.output[0]
    output:
        "{sample}/binning/concoct/cluster_attribution.tsv"
    run:
        with open(input[0]) as fin, open(output[0],'w') as fout:
            for line in fin:
                fout.write(line.replace(',','\t'))

## METABAT

rule get_metabat_depth_file:
    input:
        bam = BINNING_BAM
    output:
        temp("{sample}/binning/metabat/metabat_depth.txt")
    log:
        "{sample}/binning/metabat/metabat.log"
    conda:
        "%s/metabat.yaml" % CONDAENV
    threads:
        config['threads']
    resources:
        mem = config["java_mem"]
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output} {input.bam} &> >(tee {log})
        """



rule metabat:
    input:
        depth_file = rules.get_metabat_depth_file.output,
        contigs = BINNING_CONTIGS
    output:
        "{sample}/binning/metabat/cluster_attribution.tsv",
    params:
          sensitivity = 500 if config['binning_sensitivity'] == 'sensitive' else 200,
          min_contig_len = config["binning_min_contig_length"],
          output_prefix = "{sample}/binning/bins/bin"
    benchmark:
        "logs/benchmarks/binning/metabat/{sample}.txt"
    log:
        "logs/binning/metabat/run_metabat.txt"
    conda:
        "%s/metabat.yaml" % CONDAENV
    threads:
        config["threads"]
    resources:
        mem = config["java_mem"]
    shell:
          """
          metabat2 -i {input.contigs} \
          --abdFile {input.depth_file} \
          --minContig {params.min_contig_len} \
          --numThreads {threads} \
          --maxEdges {params.sensitivity} \
          --saveCls --noBinOut\
          -o {output} \
          &> >(tee {log})
          """



#
# localrules: MAG_analyze_metabat_clusters
# rule MAG_analyze_metabat_clusters:
#     input:
#         contigs = COMBINED_CONTIGS,
#         cluster_attribution_file = "{folder}/binning/metabat/metabat_cluster_attribution.txt".format(folder=combined_contigs_folder),
#         depth_file = "{folder}/binning/metabat_depth.txt".format(folder=combined_contigs_folder)
#     output:
#         expand("{folder}/binning/metabat/{file}", folder=combined_contigs_folder,
#                file=['cluster_attribution.txt',
#                      'contig_stats.tsv',
#                      'cluster_stats.tsv',
#                      'average_cluster_abundance.tsv',
#                      'average_contig_abundance.tsv.gz'])
#         # {folder}/binning/bins/MAG{id}.fasta
#     params:
#         output_prefix = lambda wc, output: os.path.join(os.path.dirname(output[0]), 'bins', 'Bin')
#     log:
#         "logs/binning/analyze_metabat_clusters.txt"
#     shell:
#         """
#             python %s/rules/analyze_metabat_clusters.py \
#             {input.contigs} \
#             {input.cluster_attribution_file} \
#             {input.depth_file} \
#             {params.output_prefix} \
#             2> >(tee {log})
#         """ % os.path.dirname(os.path.abspath(workflow.snakefile))
#
#
# # https://bitbucket.org/berkeleylab/metabat/wiki/Best%20Binning%20Practices
#


rule maxbin:
    input:
        fasta = BINNING_CONTIGS,
        covarage = rules.get_contig_covarage_from_bb.output
    output:
        # fastas will need to be dynamic if we do something with them at a later time
        summary = "{sample}/binning/maxbin/bins/bin.summary",
        marker = "{sample}/binning/maxbin/bins/bin.marker",
        cluster_attribution_file = "{sample}/binning/maxbin/cluster_attribution.tsv"
    params:
        mi = config["maxbin"]["max_iteration"],
        mcl = config["binning_min_contig_length"],
        pt = config["maxbin"]["prob_threshold"],
        output_prefix = lambda wildcards, output: os.path.splitext(output.summary)[0]
    log:
        "{sample}/logs/maxbin2.log"
    conda:
        "%s/optional_genome_binning.yaml" % CONDAENV
    threads:
        config["threads"]
    shell:
        """run_MaxBin.pl -contig {input.fasta} \
               -abund {input.covarage} \
               -out {params.output_prefix} \
               -min_contig_length {params.mcl} \
               -thread {threads} \
               -prob_threshold {params.pt} \
               -max_iteration {params.mi} > {log}


            cp {params.output_prefix}.marker {output.cluster_attribution_file} 2>> {log}
        """
