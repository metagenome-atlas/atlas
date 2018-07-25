

BINNING_BAM= "{sample}/sequence_alignment/{sample}.bam"
BINNING_CONTIGS= "{sample}/{sample}_contigs.fasta"
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
        min_length=config['concoct']["min_contig_length"],
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
            --iterations {params.niterations} &> >(tee {log}) 2>1
        """
localrules: concoct

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
          sensitivity = 500 if config['metabat']['sensitivity'] == 'sensitive' else 200,
          min_contig_len = config['metabat']["min_contig_length"],
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

ruleorder: maxbin > get_bins

rule maxbin:
    input:
        fasta = BINNING_CONTIGS,
        covarage = rules.get_contig_covarage_from_bb.output
    output:
        directory("{sample}/binning/maxbin/bins")

    params:
        mi = config["maxbin"]["max_iteration"],
        mcl = config["maxbin"]["min_contig_length"],
        pt = config["maxbin"]["prob_threshold"],
        output_prefix = lambda wc,output: os.path.join(output[0],wc.sample)
    log:
        "{sample}/logs/binning/maxbin.log"
    conda:
        "%s/maxbin.yaml" % CONDAENV
    threads:
        config["threads"]
    shell:
        """
            mkdir {output[0]} 2> {log}
            run_MaxBin.pl -contig {input.fasta} \
              -abund {input.covarage} \
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

rule get_maxbin_cluster_attribution:
    input:
        directory("{sample}/binning/maxbin/bins")
    output:
        "{sample}/binning/maxbin/cluster_attribution.tsv"
    params:
        file_name = lambda wc,input: "{folder}/{sample}.{{binid}}.fasta".format(folder= input[0],**wc)
    run:
        Bin_ids, = glob_wildcards(file_name)
        print("found {} bins".format(len(Bin_ids)))

        with open(output[0],'w') as out_file:

            for binid in Bin_ids:

                with open(file_name.format(binid=binid)) as bin_file:
                    for line in file:
                        if line[0] =='>':
                            fasta_header = line[1:].strip().split()[0]

                            out_file.write("{binid}/t{fasta_header}\n")




localrules: get_bins

rule get_bins:
    input:
        cluster_attribution = "{sample}/binning/{binner}/cluster_attribution.tsv",
        contigs= BINNING_CONTIGS
    output:
        directory("{sample}/binning/{binner}/bins")
    params:
        prefix= lambda wc, output: os.path.join(output[0],wc.sample)
    conda:
        "%s/sequence_utils.yaml" % CONDAENV
    script:
        "get_fasta_of_bins.py"



## Checkm

# TODO generalize checkm rules


rule initialize_checkm:
    # input:
    output:
        touched_output = "logs/checkm_init.txt"
    params:
        database_dir = CHECKMDIR
    conda:
        "%s/checkm.yaml" % CONDAENV
    log:
        "logs/initialize_checkm.log"
    shell:
        "python %s/rules/initialize_checkm.py {params.database_dir} {output.touched_output} {log}" % os.path.dirname(os.path.abspath(workflow.snakefile))




rule run_checkm_lineage_wf:
    input:
        touched_output = "logs/checkm_init.txt",
        bins = directory("{sample}/binning/{binner}/bins") # actualy path to fastas
    output:
        "{sample}/binning/{binner}/checkm/completeness.tsv"
    params:
        output_dir = lambda wc, output: os.path.dirname(output[0])
    conda:
        "%s/checkm.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """rm -r {params.output_dir} && \
           checkm lineage_wf \
               --file {params.output_dir}/completeness.tsv \
               --tab_table \
               --quiet \
               --extension fasta \
               --threads {threads} \
               {input.bins} \
               {params.output_dir}"""


rule run_checkm_tree_qa:
    input:
        "{sample}/binning/{binner}/checkm/completeness.tsv"
    output:
        "{sample}/binning/{binner}/checkm/taxonomy.tsv"
    params:
        output_dir = lambda wc,output: os.path.dirname(output[0])
    conda:
        "%s/checkm.yaml" % CONDAENV
    shell:
        """checkm tree_qa \
               --tab_table \
               --out_format 2 \
               --file {params.output_dir}/taxonomy.tsv \
               {params.output_dir}"""



rule build_bin_report:
    input:
        completeness_files = expand("{sample}/binning/{{binner}}/checkm/completeness.tsv", sample=SAMPLES),
        taxonomy_files = expand("{sample}/binning/{{binner}}/checkm/taxonomy.tsv", sample=SAMPLES)
    output:
        report = "reports/bin_report_{binner}.html",
        bin_table = "reports/genomic_bins_{binner}.tsv"
    params:
        samples = " ".join(SAMPLES)
    conda:
        "%s/report.yaml" % CONDAENV
    shell:
        """
        python %s/report/bin_report.py \
            --samples {params.samples} \
            --completeness {input.completeness_files} \
            --taxonomy {input.taxonomy_files} \
            --report-out {output.report} \
            --bin-table {output.bin_table}
        """ % os.path.dirname(os.path.abspath(workflow.snakefile))

# not working correctly https://github.com/cmks/DAS_Tool/issues/13
rule run_das_tool:
    input:
        cluster_attribution= expand("{{sample}}/binning/{binner}/cluster_attribution.tsv",
               binner=config['binner']),
        contigs = BINNING_CONTIGS,
        proteins= "{sample}/annotation/predicted_genes/{sample}.faa"
    output:
        expand("{{sample}}/binning/DASTool/{{sample}}{postfix}",
               postfix= [ "_summary.txt","_hqBins.pdf", "_scores.pdf" ]),
        cluster_attribution= "{sample}/binning/DASTool/cluster_attribution.tsv"
#        [method].eval
    threads:
        config['threads']
    log:
        "{sample}/logs/binning/DASTool.log"
    conda:
        "%s/DASTool.yaml" % CONDAENV
    params:
        binner_names= ",".join(config['binner']),
        scaffolds2bin= lambda wc, input: ",".join(input.cluster_attribution),
        output_prefix= "{sample}/binning/DASTool/{sample}",
        score_threshold = config['DASTool']['score_threshold'],
        megabin_penalty = config['DASTool']['megabin_penalty'],
        duplicate_penalty = config['DASTool']['duplicate_penalty']

    shell:
        " DAS_Tool --outputbasename {params.output_prefix} "
        " --bins {params.scaffolds2bin} "
        " --labels {params.binner_names} "
        " --contigs {input.contigs} "
        " --search_engine diamond "
        " --proteins {input.proteins} "
        " --write_bin_evals 1 "
        " --create_plots 1 --write_bin_evals 1 "
        " --megabin_penalty {params.megabin_penalty}"
        " --duplicate_penalty {params.duplicate_penalty} "
        " --threads {threads} "
        " --score_threshold {params.score_threshold} &> >(tee {log}) "
        " ; mv {params.output_prefix}_scaffolds2bin.txt {output.cluster_attribution} &> >(tee -a {log})"
