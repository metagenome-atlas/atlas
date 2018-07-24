

BINNING_BAM = "{sample}/sequence_alignment/{sample}.bam"
BINNING_CONTIGS = "{sample}/{sample}_contigs.fasta"
BB_COVERAGE_FILE = "{sample}/assembly/contig_stats/postfilter_coverage_stats.txt"


ruleorder: bam_2_sam > align_reads_to_final_contigs

rule bam_2_sam:
    input:
        "{file}.bam"
    output:
        temp("{file}.sam")
    threads:
        config['threads']
    resources:
        mem = config["java_mem"],
        java_mem = int(config["java_mem"] * JAVA_MEM_FRACTION)
    shadow:
        "shallow"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    shell:
        """
        reformat.sh in={input} out={output}
        """


localrules: get_contig_coverage_from_bb
rule get_contig_coverage_from_bb:
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
        basename = lambda wc, output: os.path.dirname(output[0]),
        Nexpected_clusters = config['concoct']['Nexpected_clusters'],
        read_length = config['concoct']['read_length'],
        min_length = config["binning_min_contig_length"],
        niterations = config["concoct"]["Niterations"]
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
        concoct -c {params.Nexpected_clusters} \
            --coverage_file {input.coverage} \
            --composition_file {input.fasta} \
            --basename {params.basename} \
            --read_length {params.read_length} \
            --length_threshold {params.min_length} \
            --converge_out \
            --iterations {params.niterations} &> >(tee {log}) 2>1
        """


localrules: convert_concoct_csv_to_tsv
rule convert_concoct_csv_to_tsv:
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
        jgi_summarize_bam_contig_depths --outputDepth {output} {input.bam} \
            &> >(tee {log})
        """


rule metabat:
    input:
        depth_file = rules.get_metabat_depth_file.output,
        contigs = BINNING_CONTIGS
    output:
        "{sample}/binning/metabat/cluster_attribution.tsv",
    params:
          sensitivity = 500 if config['metabat']['sensitivity'] == 'sensitive' else 200,
          min_contig_len = config["binning_min_contig_length"],
          output_prefix = "{sample}/binning/bins/bin"
    benchmark:
        "logs/benchmarks/binning/metabat/{sample}.txt"
    log:
        "{sample}/logs/binning/metabat.txt"
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
            --saveCls --noBinOut \
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
        covarage = rules.get_contig_coverage_from_bb.output
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
        "%s/maxbin.yaml" % CONDAENV
    threads:
        config["threads"]
    shell:
        """
        run_MaxBin.pl -contig {input.fasta} \
            -abund {input.covarage} \
            -out {params.output_prefix} \
            -min_contig_length {params.mcl} \
            -thread {threads} \
            -prob_threshold {params.pt} \
            -max_iteration {params.mi} > {log}

        cp {params.output_prefix}.marker {output.cluster_attribution_file} 2>> {log}
        """



## Checkm
# TODO generalize checkm rules
rule initialize_checkm:
    # input:
    output:
        touched_output = "logs/checkm_init.txt"
    params:
        database_dir = CHECKMDIR,
        script_dir = os.path.dirname(os.path.abspath(workflow.snakefile))
    conda:
        "%s/checkm.yaml" % CONDAENV
    log:
        "logs/initialize_checkm.log"
    shell:
        """
        python {params.script_dir}/rules/initialize_checkm.py \
            {params.database_dir} \
            {output.touched_output} \
            {log}
        """


rule run_checkm_lineage_wf:
    input:
        touched_output = "logs/checkm_init.txt",
        bins = "{sample}/binning/{binner}/bins/bin.marker"
    output:
        "{sample}/binning/checkm/{binner}/completeness.tsv"
    params:
        bin_dir = lambda wc, input: os.path.dirname(input.bins),
        output_dir = lambda wc, output: os.path.dirname(output[0])
    conda:
        "%s/checkm.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """
        rm -r {params.output_dir}
        checkm lineage_wf \
            --file {params.output_dir}/completeness.tsv \
            --tab_table \
            --quiet \
            --extension fasta \
            --threads {threads} \
            {params.bin_dir} \
            {params.output_dir}
        """


rule run_checkm_tree_qa:
    input:
        "{sample}/binning/checkm/{binner}/completeness.tsv"
    output:
        "{sample}/binning/checkm/{binner}/taxonomy.tsv"
    params:
        output_dir = "{sample}/binning/checkm/{binner}"
    conda:
        "%s/checkm.yaml" % CONDAENV
    shell:
        """
        checkm tree_qa \
            --tab_table \
            --out_format 2 \
            --file {params.output_dir}/taxonomy.tsv \
            {params.output_dir}
        """


rule build_bin_report:
    input:
        completeness_files = expand("{sample}/binning/checkm/{{binner}}/completeness.tsv", sample=SAMPLES),
        taxonomy_files = expand("{sample}/binning/checkm/{{binner}}/taxonomy.tsv", sample=SAMPLES)
    output:
        report = "reports/bin_report_{binner}.html",
        bin_table = "genomic_bins_{binner}.tsv"
    params:
        samples = " ".join(SAMPLES),
        script_dir = os.path.dirname(os.path.abspath(workflow.snakefile))
    conda:
        "%s/report.yaml" % CONDAENV
    shell:
        """
        python {params.script_dir}/report/bin_report.py \
            --samples {params.samples} \
            --completeness {input.completeness_files} \
            --taxonomy {input.taxonomy_files} \
            --report-out {output.report} \
            --bin-table {output.bin_table}
        """
