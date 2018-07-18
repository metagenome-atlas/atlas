
BINNING_BAMS="contigs/sequence_alignment_combined_contigs/{sample}.bam"
BINNING_CONTIGS
SAMPLE

config['binning_sensitivity'] = 'sensitive'
config["metabat_min_contig_length"] =500
## CONCOCT



rule run_concoct:
    input:
        coverage = "{folder}/sequence_alignment_{Reference}/combined_median_coverage.tsv".format(Reference='combined_contigs', folder=combined_contigs_folder),
        fasta = "{folder}/{Reference}.fasta".format(Reference='combined_contigs', folder=combined_contigs_folder)
    output:
        expand("{folder}/binning/{file}",
          folder=combined_contigs_folder,
            file=['means_gt2500.csv',
                'PCA_components_data_gt2500.csv',
                'original_data_gt2500.csv',
                'PCA_transformed_data_gt2500.csv',
                'pca_means_gt2500.csv',
                'args.txt',
                'responsibilities.csv']
        ),
    params:
        basename= lambda wc,output: os.path.dirname(output[0]),
          Nexpected_clusters= config['concoct']['Nexpected_clusters'],
          read_length= config['concoct']['read_length'],
          min_length=config["minimum_contig_length"],
          niterations=config["concoct"]["Niterations"]
    benchmark:
        "logs/benchmarks/binning/concoct.txt"
    log:
        "{folder}/binning/log.txt".format(folder=combined_contigs_folder)
    conda:
        "%s/concoct.yaml" % CONDAENV
    threads:
        10 # concoct uses 10 threads by default, wit for update: https://github.com/BinPro/CONCOCT/issues/177
    resources:
        mem = config.get("java_mem", JAVA_MEM)
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


rule get_metabat_depth_file:
    input:
        bam = expand(BINNING_BAMS, sample=SAMPLES)
    output:
        "{folder}/binning/metabat/metabat_depth.txt"
    log:
        "{folder}/binning/metabat/metabat.log"
    conda:
        "%s/metabat.yaml" % CONDAENV
    threads:
        config['threads']
    resources:
        mem = config.get("java_mem", JAVA_MEM)
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output} {input.bam} &> >(tee {log})
        """

    rule run_metabat:
        input:
            depth_file = "{folder}/binning/metabat/metabat_depths.txt",
            contigs = CONTIGS_FOR_BINNING
        output:
            "{folder}/binning/metabat/metabat_cluster_attribution.txt",
        params:
              sensitivity = 500 if config['binning_sensitivity'] == 'sensitive' else 200,
              min_contig_len = config["metabat_min_contig_length"],
              output_prefix = "{folder}/binning/bins/bin"
        benchmark:
            "logs/benchmarks/binning/metabat.txt"
        log:
            "logs/binning/metabat/run_metabat.txt"
        conda:
            "%s/metabat.yaml" % CONDAENV
        threads:
            config["threads"]
        resources:
            mem = config.get("java_mem", JAVA_MEM)
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

    localrules: MAG_analyze_metabat_clusters
    rule MAG_analyze_metabat_clusters:
        input:
            contigs = COMBINED_CONTIGS,
            cluster_attribution_file = "{folder}/binning/metabat/metabat_cluster_attribution.txt".format(folder=combined_contigs_folder),
            depth_file = "{folder}/binning/metabat_depth.txt".format(folder=combined_contigs_folder)
        output:
            expand("{folder}/binning/metabat/{file}", folder=combined_contigs_folder,
                   file=['cluster_attribution.txt',
                         'contig_stats.tsv',
                         'cluster_stats.tsv',
                         'average_cluster_abundance.tsv',
                         'average_contig_abundance.tsv.gz'])
            # {folder}/binning/bins/MAG{id}.fasta
        params:
            output_prefix = lambda wc, output: os.path.join(os.path.dirname(output[0]), 'bins', 'Bin')
        log:
            "logs/binning/analyze_metabat_clusters.txt"
        shell:
            """
                python %s/rules/analyze_metabat_clusters.py \
                {input.contigs} \
                {input.cluster_attribution_file} \
                {input.depth_file} \
                {params.output_prefix} \
                2> >(tee {log})
            """ % os.path.dirname(os.path.abspath(workflow.snakefile))


# https://bitbucket.org/berkeleylab/metabat/wiki/Best%20Binning%20Practices


    localrules: make_maxbin_abundance_list
    rule make_maxbin_abundance_list:
        input:
            coverage = expand("contigs/sequence_alignment_combined_contigs/{sample}/{sample}_coverage.txt",
                sample=SAMPLES)
        output:
            abundance_file = expand("contigs/maxbin/contig_coverage/{sample}.tsv", sample=SAMPLES),
            abundance_list = temp("contigs/maxbin/contig_coverage.list")
        run:
            with open(output.abundance_list, 'w') as outf:
                for i in range(len(input)):
                    bb_cov_stats_to_maxbin(input.coverage[i], output.abundance_file[i])
                    outf.write(os.path.abspath(output.abundance_file[i]) + '\n')


    rule MAG_run_maxbin:
        input:
            fasta = COMBINED_CONTIGS,
            abundance_list = "contigs/maxbin/contig_coverage.list"
        output:
            # fastas will need to be dynamic if we do something with them at a later time
            summary = "contigs/maxbin/bins/bin.summary",
            marker = "contigs/maxbin/bins/bin.marker",
            cluster_attribution_file = "contigs/maxbin/cluster_attribution.txt"
        benchmark:
            "logs/benchmarks/maxbin2/combined_contigs.txt"
        params:
            mi = config.get("maxbin_max_iteration", MAXBIN_MAX_ITERATION),
            mcl = config.get("maxbin_min_contig_length", MAXBIN_MIN_CONTIG_LENGTH),
            pt = config.get("maxbin_prob_threshold", MAXBIN_PROB_THRESHOLD),
            out_prefix = lambda wildcards, output: os.path.splitext(os.path.dirname(output.summary))[0]
        log:
            "contigs/logs/maxbin2.log"
        conda:
            "%s/optional_genome_binning.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        shell:
            """run_MaxBin.pl -contig {input.fasta} \
                   -abund_list {input.abundance_list} \
                   -out {params.outdir} \
                   -min_contig_length {params.mcl} \
                   -thread {threads} \
                   -prob_threshold {params.pt} \
                   -max_iteration {params.mi} > {log}


                cp {params.output_prefix}.marker {input.cluster_attribution_file} 2>> {log}
            """






    # rule MAG_initialize_checkm:
    #     # input:
    #     output:
    #         touched_output = "logs/checkm_init.txt"
    #     params:
    #         database_dir = CHECKMDIR
    #     conda:
    #         "%s/optional_genome_binning.yaml" % CONDAENV
    #     log:
    #         "logs/initialize_checkm.log"
    #     shell:
    #         "python %s/rules/initialize_checkm.py {params.database_dir} {output.touched_output} {log}" % os.path.dirname(os.path.abspath(workflow.snakefile))


    rule MAG_run_checkm_lineage_wf:
        input:
            touched_output = "logs/checkm_init.txt",
            bins = "{folder}/binning/cluster_attribution.txt".format(folder=combined_contigs_folder)
        output:
            "{folder}/binning/checkm/completeness.tsv".format(folder=combined_contigs_folder)
        params:
            bin_dir = lambda wc, input: os.path.join(os.path.dirname(input.bins),"bins"),
            output_dir = lambda wc, output: os.path.dirname(output[0]),
            fasta_extension = 'fasta'
        conda:
            "%s/optional_genome_binning.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        resources:
            mem= config.get('java_mem',JAVA_MEM)
        shell:
            """rm -r {params.output_dir} && \
               checkm lineage_wf \
                   --file {params.output_dir}/completeness.tsv \
                   --tab_table \
                   --quiet \
                   --extension {params.fasta_extension} \
                   --threads {threads} \
                   {params.bin_dir} \
                   {params.output_dir}"""


    rule MAG_run_checkm_tree_qa:
        input:
            "{folder}/binning/checkm/completeness.tsv".format(folder=combined_contigs_folder)
        output:
            "{folder}/binning/checkm/taxonomy.tsv".format(folder=combined_contigs_folder)
        params:
            output_dir = lambda wc, output: os.path.dirname(output[0])
        conda:
            "%s/optional_genome_binning.yaml" % CONDAENV
        shell:
            """checkm tree_qa \
                   --tab_table \
                   --out_format 2 \
                   --file {params.output_dir}/taxonomy.tsv \
                   {params.output_dir}"""
