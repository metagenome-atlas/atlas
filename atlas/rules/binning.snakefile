
from glob import glob
BINNING_CONTIGS= "{sample}/{sample}_contigs.fasta"


rule bam_2_sam_binning:
    input:
        "{sample}/sequence_alignment/{sample_reads}.bam"
    output:
        temp("{sample}/sequence_alignment/binning_{sample_reads}.sam")
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
        reformat.sh in={input} out={output} sam=1.3
        """




rule pileup_for_binning:
    input:
        fasta = BINNING_CONTIGS,
        sam = rules.bam_2_sam_binning.output,
    output:
        covstats = "{sample}/binning/coverage/{sample_reads}_coverage_stats.txt",
    params:
        pileup_secondary = 't' if config.get("count_multi_mapped_reads", CONTIG_COUNT_MULTI_MAPPED_READS) else 'f',
    log:
        "{sample}/logs/binning/calculate_coverage/pileup_reads_from_{sample_reads}_to_filtered_contigs.log" # this file is udes for assembly report
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", JAVA_MEM),
        java_mem = int(config.get("java_mem", JAVA_MEM) * JAVA_MEM_FRACTION)
    shell:
        """pileup.sh ref={input.fasta} in={input.sam} \
               threads={threads} \
               -Xmx{resources.java_mem}G \
               covstats={output.covstats} \
               secondary={params.pileup_secondary} \
                2> {log}
        """




localrules: get_contig_coverage_from_bb, combine_coverages
rule get_contig_coverage_from_bb:
    input:
        coverage = "{sample}/binning/coverage/{sample_reads}_coverage_stats.txt"
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
        covstats = lambda wc: expand("{sample}/binning/coverage/{sample_reads}_coverage_stats.txt",
                                 sample_reads = GROUPS[config['samples'][wc.sample]['group']],
                                 sample=wc.sample)
    output:
        "{sample}/binning/coverage/combined_coverage.tsv"
    run:

        import pandas as pd
        import os

        combined_cov={}
        for cov_file in input:

            sample= os.path.split(cov_file)[-1].split('_')[0]
            data= pd.read_table(cov_file,index_col=0)

            data.loc[data.Avg_fold<0,'Avg_fold']=0
            combined_cov[sample]= data.Avg_fold
        pd.DataFrame(combined_cov).to_csv(output[0],sep='\t')



## CONCOCT
rule run_concoct:
    input:
        coverage = "{sample}/binning/coverage/combined_coverage.tsv",
        fasta = BINNING_CONTIGS
    output:
        "{{sample}}/binning/concoct/intermediate_files/clustering_gt{}.csv".format(config['concoct']["min_contig_length"])
    params:
        basename= lambda wc, output: os.path.dirname(output[0]),
        Nexpected_clusters= config['concoct']['Nexpected_clusters'],
        read_length= config['concoct']['read_length'],
        min_length=config['concoct']["min_contig_length"],
        niterations=config["concoct"]["Niterations"]
    log:
        "{sample}/binning/concoct/intermediate_files/log.txt"
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
            --iterations {params.niterations}
        """


localrules: convert_concoct_csv_to_tsv
rule convert_concoct_csv_to_tsv:
    input:
        rules.run_concoct.output[0]
    output:
        "{sample}/binning/concoct/cluster_attribution.tmp"
    run:
        with open(input[0]) as fin, open(output[0],'w') as fout:
            for line in fin:
                fout.write(line.replace(',','\t'))


## METABAT
rule get_metabat_depth_file:
    input:
        bam = lambda wc: expand("{sample}/sequence_alignment/{sample_reads}.bam",
                     sample_reads = GROUPS[config['samples'][wc.sample]['group']],
                     sample=wc.sample)
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
        "{sample}/binning/metabat/cluster_attribution.tmp",
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
            --saveCls --noBinOut \
            -o {output} \
            &> >(tee {log})
        """





rule maxbin:
    input:
        fasta = BINNING_CONTIGS,
        abund = "{sample}/binning/coverage/{sample}_coverage.txt",
    output:
        directory("{sample}/binning/maxbin/intermediate_files")
    params:
        mi = config["maxbin"]["max_iteration"],
        mcl = config["maxbin"]["min_contig_length"],
        pt = config["maxbin"]["prob_threshold"],
        output_prefix = lambda wc, output: os.path.join(output[0], wc.sample)
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





localrules: get_maxbin_cluster_attribution, get_bins
localrules: get_unique_cluster_attribution,
rule get_unique_cluster_attribution:
    input:
        "{sample}/binning/{binner}/cluster_attribution.tmp"
    output:
        "{sample}/binning/{binner}/cluster_attribution.tsv"
    run:
        import pandas as pd
        import numpy as np


        d= pd.read_table(input[0],index_col=0, squeeze=True, header=None)

        assert type(d) == pd.Series, "expect the input to be a two column file: {}".format(input[0])

        old_cluster_ids = list(d.unique())
        if 0 in old_cluster_ids:
            old_cluster_ids.remove(0)

        map_cluster_ids = dict(zip(old_cluster_ids, gen_names_for_range(len(old_cluster_ids), prefix="{sample}_{binner}_" )  ))

        new_d= d.map(map_cluster_ids)
        new_d.dropna(inplace=True)

        new_d.to_csv(output[0],sep='\t')


rule get_maxbin_cluster_attribution:
    input:
        directory("{sample}/binning/maxbin/intermediate_files")
    output:
        "{sample}/binning/maxbin/cluster_attribution.tsv"
    params:
        file_name = lambda wc, input: "{folder}/{sample}.{{binid}}.fasta".format(folder=input[0], **wc)
    run:
        bin_ids, = glob_wildcards(params.file_name)
        print("found {} bins".format(len(bin_ids)))
        with open(output[0],'w') as out_file:
            for binid in bin_ids:
                with open(params.file_name.format(binid=binid)) as bin_file:
                    for line in bin_file:
                        if line.startswith(">"):
                            fasta_header = line[1:].strip().split()[0]
                            out_file.write("{fasta_header}\t{sample}_maxbin_{binid}\n".format(binid=binid,
                                                                                              fasta_header=fasta_header,
                                                                                              sample=wildcards.sample))
                os.remove(params.file_name.format(binid=binid))


rule get_bins:
    input:
        cluster_attribution = "{sample}/binning/{binner}/cluster_attribution.tsv",
        contigs= BINNING_CONTIGS
    output:
        directory("{sample}/binning/{binner}/bins")
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
        """
        rm -r {params.output_dir}
        checkm lineage_wf \
            --file {params.output_dir}/completeness.tsv \
            --tab_table \
            --quiet \
            --extension fasta \
            --threads {threads} \
            {input.bins} \
            {params.output_dir}
        """



rule run_checkm_tree_qa:
    input:
        tree="{checkmfolder}/completeness.tsv"
    output:
        netwick="{checkmfolder}/tree.nwk",
        summary="{checkmfolder}/taxonomy.tsv",
    params:
        tree_dir = lambda wc, input: os.path.dirname(input.tree),
    conda:
        "%s/checkm.yaml"  % CONDAENV
    threads:
        1
    shell:
        """
            checkm tree_qa \
               {params.tree_dir} \
               --out_format 4 \
               --file {output.netwick}

               checkm tree_qa \
                  {params.tree_dir} \
                  --out_format 2 \
                  --file {output.summary}\
                  --tab_table

        """


rule checkm_tetra:
    input:
        contigs=BINNING_CONTIGS,
    output:
        "{sample}/binning/{binner}/checkm/tetranucleotides.txt"
    log:
        "{sample}/logs/binning/{binner}/checkm/tetra.txt"
    conda:
        "%s/checkm.yaml" % CONDAENV
    threads:
        config.get("threads", 8)
    shell:
        """
            checkm tetra \
            --threads {threads} \
            {input.contigs} {output} 2> {log}
        """


rule checkm_outliers:
    input:
        tetra= "{sample}/binning/{binner}/checkm/tetranucleotides.txt",
        bin_folder= directory("{sample}/binning/{binner}/bins"),
        checkm = "{sample}/binning/{binner}/checkm/completeness.tsv"
    params:
        checkm_folder = lambda wc, input: os.path.dirname(input.checkm),
        report_type = 'any',
        treshold = 95 #reference distribution used to identify outliers; integer between 0 and 100 (default: 95)
    output:
        "{sample}/binning/{binner}/checkm/outliers.txt"
    log:
        "{sample}/logs/binning/{binner}/checkm/outliers.txt"
    conda:
        "%s/checkm.yaml" % CONDAENV
    threads:
        config.get("threads", 8)
    shell:
        """
            checkm outliers \
            --extension fasta \
            --distributions {params.treshold} \
            --report_type {params.report_type} \
            {params.checkm_folder} \
            {input.bin_folder} \
            {input.tetra} \
            {output} 2> {log}
        """


rule find_16S:
    input:
        contigs=BINNING_CONTIGS,
        bin_dir= directory("{sample}/binning/{binner}/bins")
    output:
        summary="{sample}/binning/{binner}/SSU/ssu_summary.tsv",
        fasta="{sample}/binning/{binner}/SSU/ssu.fna",
    params:
        output_dir = lambda wc, output: os.path.dirname(output[0]),
        evalue = 1e-05,
        concatenate = 200 #concatenate hits that are within the specified number of base pairs
    conda:
        "%s/checkm.yaml" % CONDAENV
    threads:
        1
    shell:
        """
        rm -r {params.output_dir} && \
           checkm ssu_finder \
               --extension fasta \
               --threads {threads} \
               --evalue {params.evalue} \
               --concatenate {params.concatenate} \
               {input.contigs} \
               {input.bin_dir} \
               {params.output_dir}
        """


rule get_all_16S:
    input:
        summaries= expand(rules.find_16S.output.summary,sample=SAMPLES,binner=config['final_binner']),
        fastas= expand(rules.find_16S.output.fasta,sample=SAMPLES,binner=config['final_binner'])
    output:
        fasta="genomes/SSU/ssu.fasta",
        summary ="genomes/SSU/ssu_summary.tsv"
    run:
        shell("cat {input.fastas} > {output.fasta}")

        import pandas as pd
        summary= pd.DataFrame()

        for file in input.summaries:
            try:
                d = pd.read_table(file,index_col=0)
                summary=summary.append(d)
            except:
                pd.errors.EmptyDataError

        summary.to_csv(output.summary,sep='\t')




localrules: build_bin_report
rule build_bin_report:
    input:
        completeness_files = expand("{sample}/binning/{{binner}}/checkm/completeness.tsv", sample=SAMPLES),
        taxonomy_files = expand("{sample}/binning/{{binner}}/checkm/taxonomy.tsv", sample=SAMPLES)
    output:
        report = "reports/bin_report_{binner}.html",
        bin_table = "reports/genomic_bins_{binner}.tsv"
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

localrules: get_unique_bin_ids
rule get_unique_bin_ids:
    input:
        "{sample}/binning/{binner}/cluster_attribution.tsv"
    output:
        "{sample}/binning/DASTool/{binner}.scaffolds2bin"
    shell:
        "cp {input} {output}"


rule run_das_tool:
    input:
        cluster_attribution = expand("{{sample}}/binning/DASTool/{binner}.scaffolds2bin",
            binner=config['binner']),
        contigs = BINNING_CONTIGS,
        proteins= "{sample}/annotation/predicted_genes/{sample}.faa"
    output:
        expand("{{sample}}/binning/DASTool/{{sample}}{postfix}",
               postfix=["_DASTool_summary.txt", "_DASTool_hqBins.pdf", "_DASTool_scores.pdf"]),
        expand("{{sample}}/binning/DASTool/{{sample}}_{binner}.eval",
               binner= config['binner']),
        cluster_attribution = "{sample}/binning/DASTool/cluster_attribution.tsv"
    threads:
        config['threads']
    log:
        "{sample}/logs/binning/DASTool.log"
    conda:
        "%s/DASTool.yaml" % CONDAENV
    params:
        binner_names = ",".join(config['binner']),
        scaffolds2bin = lambda wc, input: ",".join(input.cluster_attribution),
        output_prefix = "{sample}/binning/DASTool/{sample}",
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
        " --debug "
        " --score_threshold {params.score_threshold} &> >(tee {log}) "
        " ; mv {params.output_prefix}_DASTool_scaffolds2bin.txt {output.cluster_attribution} &> >(tee -a {log})"


# # unknown bins and contigs
#
if config['final_binner']=='DASTool':
    localrules: get_unknown_bins
    rule get_unknown_bins:
        input:
            score_files=expand("{{sample}}/binning/DASTool/{{sample}}_{binner}.eval", binner= config['binner']),
            bin_dirs=expand(directory("{{sample}}/binning/{binner}/bins"), binner= config['binner']),
        output:
            dir= directory("{sample}/binning/Unknown/bins"),
            scores= "{sample}/binning/Unknown/scores.tsv"

        run:
            import pandas as pd
            import shutil

            Scores= pd.DataFrame()
            for i in range(len(config['binner'])):
                score_file = input.score_files[i]
                bin_dir = input.bin_dirs[i]

                S = pd.read_table(score_file,index_col=0)

                S= S.loc[S.SCG_completeness==0]
                Scores= Scores.append(S)

                for bin_id in S.index:
                    shutil.copy(os.path.join(bin_dir,bin_id+'.fasta'), output.dir )

            Scores.to_csv(output.scores,sep='\t')
    #
    #
    #


## dRep
localrules: get_all_bins
rule get_all_bins:
    input:
        bins=expand(directory("{sample}/binning/{binner}/bins"),
               sample= SAMPLES, binner= config['final_binner']),
        cluster_attribution=expand("{sample}/binning/{binner}/cluster_attribution.tsv",
               sample= SAMPLES, binner= config['final_binner'])
    output:
        temp(directory("genomes/all_bins"))

    run:
        os.mkdir(output[0])
        from glob import glob
        import shutil
        for bin_folder in input.bins:
            for fasta_file in glob(bin_folder+'/*.fasta'):

                #fasta_file_name = os.path.split(fasta_file)[-1]
                #in_path = os.path.dirname(fasta_file)
                #out_path= os.path.join(output[0],fasta_file_name)
                #os.symlink(os.path.relpath(fasta_file,output[0]),out_path)

                shutil.copy(fasta_file,output[0])



localrules: get_quality_for_dRep_from_checkm
rule get_quality_for_dRep_from_checkm:
    input:
        "reports/genomic_bins_{binner}.tsv".format(binner=config['final_binner'])
    output:
        "genomes/quality.csv"
    run:
        import pandas as pd

        D= pd.read_table(input[0],index_col=0)

        D.index+=".fasta"
        D.index.name="genome"
        D.columns= D.columns.str.lower()
        D.iloc[:,:3].to_csv(output[0])

if config['final_binner']=='DASTool':
    localrules: get_quality_for_dRep_from_DASTool
    ruleorder: get_quality_for_dRep_from_DASTool> get_quality_for_dRep_from_checkm
    rule get_quality_for_dRep_from_DASTool:
        input:
            expand("{sample}/binning/DASTool/{sample}_DASTool_summary.txt",sample=SAMPLES)
        output:
            "genomes/DASTool_quality.tsv",
            "genomes/quality.csv"
        run:
            import pandas as pd
            D= pd.DataFrame()
            for i,file in enumerate(input):
                d= pd.read_table(file,index_col=0)
                d.index= SAMPLES[i]+'.'+d.index
                D= D.append(d)

            D.to_csv(output[0],sep='\t')

            D.index+=".fasta"
            D.index.name="genome"

            D= D.rename(columns={"SCG_completeness":"completeness", "SCG_redundancy":"contamination"})

            D[["completeness","contamination"]].to_csv(output[1])

rule first_dereplication:
    input:
        directory("genomes/all_bins"),
        quality= "genomes/quality.csv"
    output:
        directory("genomes/pre_dereplication/dereplicated_genomes")
    threads:
        config['threads']
    log:
        "logs/genomes/pre_dereplication.log"
    conda:
        "%s/dRep.yaml" % CONDAENV
    params:
        filter= " --noQualityFiltering " if config['genome_dereplication']['filter']['noFilter'] else "",
        filter_length= config['genome_dereplication']['filter']['length'],
        filter_completeness= config['genome_dereplication']['filter']['completeness'],
        filter_contamination= config['genome_dereplication']['filter']['contamination'],
        ANI= config['genome_dereplication']['ANI'],
        completeness_weight= config['genome_dereplication']['score']['completeness'] ,
        contamination_weight=config['genome_dereplication']['score']['contamination'] ,
        strain_heterogeneity_weight= config['genome_dereplication']['score']['completeness'] , #not in table
        N50_weight=config['genome_dereplication']['score']['N50'] ,
        size_weight=config['genome_dereplication']['score']['length'] ,
        opt_parameters = config['genome_dereplication']['opt_parameters'],
        work_directory= lambda wc,output: os.path.dirname(output[0]),
        sketch_size= config['genome_dereplication']['sketch_size']

    shell:
        " dRep dereplicate "
        " {params.filter} "
        " --genomes {input[0]}/*.fasta "
        " --genomeInfo {input.quality} "
        " --length {params.filter_length} "
        " --completeness {params.filter_completeness} "
        " --contamination {params.filter_contamination} "
        " --SkipSecondary "
        " --P_ani {params.ANI} "
        " --completeness_weight {params.completeness_weight} "
        " --contamination_weight {params.contamination_weight} "
        " --strain_heterogeneity_weight {params.strain_heterogeneity_weight} "
        " --N50_weight {params.N50_weight} "
        " --size_weight {params.size_weight} "
        " --MASH_sketch {params.sketch_size} "
        " --processors {threads} "
        " {params.opt_parameters} "
        " {params.work_directory} "
        " &> {log} "

rule second_dereplication:
    input:
        rules.first_dereplication.output,
        quality= "genomes/quality.csv"
    output:
        directory("genomes/Dereplication/dereplicated_genomes")
    threads:
        config['threads']
    log:
        "logs/genomes/dereplication.log"
    conda:
        "%s/dRep.yaml" % CONDAENV
    params:
        ANI= config['genome_dereplication']['ANI'],
        completeness_weight= config['genome_dereplication']['score']['completeness'] ,
        contamination_weight=config['genome_dereplication']['score']['contamination'] ,
        strain_heterogeneity_weight= config['genome_dereplication']['score']['completeness'] , #not in table
        N50_weight=config['genome_dereplication']['score']['N50'] ,
        size_weight=config['genome_dereplication']['score']['length'] ,
        opt_parameters = config['genome_dereplication']['opt_parameters'],
        work_directory= lambda wc,output: os.path.dirname(output[0]),
        sketch_size= config['genome_dereplication']['sketch_size']

    shell:
        " dRep dereplicate "
        " --genomes {input[0]}/*.fasta "
        " --genomeInfo {input.quality} "
        " --noQualityFiltering "
        " --S_ani {params.ANI} "
        " --completeness_weight {params.completeness_weight} "
        " --contamination_weight {params.contamination_weight} "
        " --strain_heterogeneity_weight {params.strain_heterogeneity_weight} "
        " --N50_weight {params.N50_weight} "
        " --size_weight {params.size_weight} "
        " --MASH_sketch {params.sketch_size} "
        " --processors {threads} "
        " {params.opt_parameters} "
        " {params.work_directory} "
        " &> {log} "



rule run_all_checkm_lineage_wf:
    input:
        touched_output = "logs/checkm_init.txt",
        bins = directory("genomes/genomes")
    output:
        "genomes/checkm/completeness.tsv"
    params:
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
            {input.bins} \
            {params.output_dir}
        """

localrules: rename_final_bins
rule rename_final_bins:
    input:
        directory("genomes/Dereplication/dereplicated_genomes")
    output:
        dir=directory("genomes/genomes"),
        map_file="genomes/contig2genome.tsv",
        old2new = "genomes/old2newID.tsv"
    params:
        file_name = lambda wc, input: "{folder}/{{binid}}.fasta".format(folder=input[0], **wc)
    run:
        import shutil
        bin_ids, = glob_wildcards(params.file_name)

        old2new_name= dict(zip(bin_ids,gen_names_for_range(len(bin_ids),prefix='MAG')))
        os.makedirs(output.dir)

        with open(output.map_file,'w') as out_file, open(output.old2new,'w') as old2new_mapping_file :
            old2new_mapping_file.write(f"BinID\tMAG\n")
            for binid in bin_ids:

                fasta_in = params.file_name.format(binid=binid)
                new_name= old2new_name[binid]

                old2new_mapping_file.write(f"{binid}\t{new_name}\n")

                fasta_out = os.path.join(output.dir,f"{new_name}.fasta")
                shutil.copy(fasta_in,fasta_out)

                # write names of contigs in mapping file
                with open(fasta_in) as bin_file:
                    for line in bin_file:
                        if line[0]==">":
                            contig = line[1:].strip().split()[0]
                            out_file.write(f"{contig}\t{new_name}\n")


rule build_db_genomes:
    input:
        fasta_dir = directory("genomes/genomes")
    output:
        index="ref/genome/3/summary.txt",
        fasta=temp("genomes/all_contigs.fasta")
    threads:
        config.get("threads", 6)
    resources:
        mem = config["java_mem"],
        java_mem = int(config["java_mem"] * JAVA_MEM_FRACTION)
    log:
        "logs/genomes/mapping/build_bbmap_index.log"
    shell:
        """
        cat {input.fasta_dir}/* > {output.fasta} 2> {log}
        bbmap.sh build=3 -Xmx{resources.java_mem}G ref={output.fasta} threads={threads} local=f 2>> {log}

        """





# generalized rule so that reads from any "sample" can be aligned to contigs from "sample_contigs"
rule align_reads_to_MAGs:
    input:
        unpack(get_quality_controlled_reads),
        ref = rules.build_db_genomes.output.index,
    output:
        sam = temp("genomes/alignments/{sample}.sam"),
    params:
        input = lambda wc, input : input_params_for_bbwrap(wc, input),
        maxsites = config.get("maximum_counted_map_sites", MAXIMUM_COUNTED_MAP_SITES),
        max_distance_between_pairs = config.get('contig_max_distance_between_pairs', CONTIG_MAX_DISTANCE_BETWEEN_PAIRS),
        paired_only = 't' if config.get("contig_map_paired_only", CONTIG_MAP_PAIRED_ONLY) else 'f',
        ambiguous = 'all' if CONTIG_COUNT_MULTI_MAPPED_READS else 'best',
        min_id = config.get('contig_min_id', CONTIG_MIN_ID),
        maxindel = 100 # default 16000 good for genome deletions but not necessarily for alignment to contigs
    shadow:
        "shallow"
    log:
        "logs/genomes/mapping/map_{sample}.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", JAVA_MEM),
        java_mem = int(config.get("java_mem", JAVA_MEM) * JAVA_MEM_FRACTION)
    shell:
        """
            bbwrap.sh \
            build=3 \
            {params.input} \
            trimreaddescriptions=t \
            out={output.sam} \
            threads={threads} \
            pairlen={params.max_distance_between_pairs} \
            pairedonly={params.paired_only} \
            minid={params.min_id} \
            mdtag=t \
            xstag=fs \
            nmtag=t \
            sam=1.3 \
            ambiguous={params.ambiguous} \
            secondary=t \
            saa=f \
            maxsites={params.maxsites} \
            -Xmx{resources.java_mem}G \
            2> {log}
        """

ruleorder: bam_2_sam_MAGs > align_reads_to_MAGs
rule bam_2_sam_MAGs:
    input:
        "genomes/alignments/{sample}.bam"
    output:
        temp("genomes/alignments/{sample}.sam")
    threads:
        config['threads']
    resources:
        mem = config["java_mem"],
    shadow:
        "shallow"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    shell:
        """
        reformat.sh in={input} out={output} sam=1.3
        """



rule pileup_MAGs:
    input:
        sam = "genomes/alignments/{sample}.sam",
        #bam = "genomes/alignments/{sample}.bam" # to store it
    output:
        basecov = temp("genomes/alignments/{sample}_base_coverage.txt.gz"),
        covhist = temp("genomes/alignments/{sample}_coverage_histogram.txt"),
        covstats = temp("genomes/alignments/{sample}_coverage.txt"),
        bincov = temp("genomes/alignments/{sample}_coverage_binned.txt")
    log:
        "logs/genomes/alignments/pilup_{sample}.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", JAVA_MEM),
        java_mem = int(config.get("java_mem", JAVA_MEM) * JAVA_MEM_FRACTION)
    shell:
        """pileup.sh in={input.sam} \
               threads={threads} \
               -Xmx{resources.java_mem}G \
               covstats={output.covstats} \
               hist={output.covhist} \
               basecov={output.basecov}\
               concise=t \
               bincov={output.bincov} 2> {log}"""



localrules: combine_coverages_MAGs,combine_bined_coverages_MAGs
rule combine_coverages_MAGs:
    input:
        covstats = expand("genomes/alignments/{sample}_coverage.txt",
            sample=SAMPLES)
    output:
        "genomes/counts/median_contig_coverage.tsv",
        "genomes/counts/raw_counts_contigs.tsv",
    run:

        import pandas as pd
        import os

        combined_cov={}
        combined_N_reads={}
        for cov_file in input:

            sample= os.path.split(cov_file)[-1].split('_')[0]
            data= pd.read_table(cov_file,index_col=0)
            data.loc[data.Median_fold<0,'Median_fold']=0
            combined_cov[sample]= data.Median_fold
            combined_N_reads[sample] = data.Plus_reads+data.Minus_reads

        pd.DataFrame(combined_cov).to_csv(output[0],sep='\t')
        pd.DataFrame(combined_N_reads).to_csv(output[1],sep='\t')



rule combine_bined_coverages_MAGs:
    input:
        binned_coverage_files = expand("genomes/alignments/{sample}_coverage_binned.txt",
            sample=SAMPLES),
        cluster_attribution = "genomes/contig2genome.tsv"
    params:
        samples= SAMPLES
    output:
        binned_cov= "genomes/counts/binned_coverage.tsv.gz",
        median_abund = "genomes/counts/median_coverage_genomes.tsv"
    run:

        import pandas as pd
        import os

        def read_coverage_binned(covarage_binned_file):
            return pd.read_table(covarage_binned_file,
                             skiprows=2,
                             index_col=[0,2],
                             usecols=[0,1,2],
                             squeeze=True)


        binCov={}
        for i, cov_file in enumerate(input.binned_coverage_files):

            sample= params.samples[i]

            binCov[sample] = read_coverage_binned(cov_file)

        binCov = pd.DataFrame(binCov)
        binCov.index.names=['Contig','Position']
        binCov.to_csv(output.binned_cov,sep='\t',compression='gzip')

        cluster_attribution = pd.read_table(input.cluster_attribution,header=None,index_col=0,squeeze=True)

        Median_abund= binCov.groupby(cluster_attribution.loc[binCov.index.get_level_values(0)].values).median().T

        Median_abund.to_csv(output.median_abund,sep='\t')
