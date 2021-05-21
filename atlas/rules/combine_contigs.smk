import os
import re
import sys
import tempfile
import warnings

from default_values import *

combined_contigs_folder='contigs'

rule combine_contigs_report:
    input:
        combined_contigs= COMBINED_CONTIGS,
        combined_contigs_stats="contigs/combined_contigs_stats.txt",
        gene_counts= 'contigs/combined_gene_counts.tsv', gene_info= 'contigs/combined_gene_info.tsv',
        # annotation= "contigs/annotations.txt",
        #median_coverage="{folder}/sequence_alignment_{Reference}/combined_median_coverage.tsv".format(Reference='combined_contigs',folder=combined_contigs_folder),
        #gc_stats = "{folder}/combined_contigs_stats_gc.tsv".format(folder=combined_contigs_folder),
        #binned_coverage = "{folder}/sequence_alignment_{Reference}/combined_coverage_binned.tsv.gz".format(Reference='combined_contigs',folder=combined_contigs_folder),
        #concoct="{folder}/binning/{file}".format(folder=combined_contigs_folder,file='means_gt2500.csv')
        #bam=expand("{folder}/sequence_alignment_{Reference}/{sample}/{sample}.bam",sample=SAMPLES,Reference='combined_contigs',folder=combined_contigs_folder),
    output:
        touch("Combined_contigs_done")

  conf["combine_contigs"] = True
  config["combine_contigs_params"]=dict(min_overlap = 200,
                                 max_substitutions=4,
                                 dont_allow_N=True,
                                 remove_cycles=True,
                                 trim_contradictions=True, #False
                                 binner='concoct'
                                 )


#### combine contigs

config['perform_genome_binning']= False

rule combine_contigs:
    input:
        expand("{sample}/{sample}_contigs.fasta",sample=SAMPLES)
    output:
        combined_contigs=temp("{folder}/combined_contigs_oldnames.fasta"),
        cluster_stats="{folder}/combined_contigs_kmerfreq.txt",
        dot="{folder}/combined_contigs_graph.dot"
    benchmark:
        "logs/benchmarks/combine_contigs.txt"
    log:
        "logs/combine_contigs.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config["mem"]
    params:
       input=lambda wc,input: ','.join(input),
       min_length=config.get("minimum_contig_length", MINIMUM_CONTIG_LENGTH),
       min_overlap=config['combine_contigs_params']['min_overlap'],
       max_substitutions=config['combine_contigs_params']['max_substitutions'],
       dont_allow_N= 't' if config['combine_contigs_params']['dont_allow_N'] else 'f',
       remove_cycles='t' if config['combine_contigs_params']['remove_cycles'] else 'f',
       trim_contradictions='t' if config['combine_contigs_params']['trim_contradictions'] else 'f'

    shell:
        """
            dedupe.sh in={params.input} findoverlaps cluster processclusters \
            out={output.combined_contigs} \
            csf={output.cluster_stats} \
            dot={output.dot} \
            minoverlap={params.min_overlap}\
            minscaf={params.min_length} \
            maxsubs={params.max_substitutions} \
            threads={threads} \
            sort=length \
            maxspanningtree={params.remove_cycles} \
            exact={params.dont_allow_N}\
            fixcanoncontradictions={params.trim_contradictions}\
            -Xmx{resources.mem}G 2> {log}
        """

# vizualize dot, takes enormous times
# dot -Tpdf combined_contigs_graph.dot -o combined_clusters.pdf

localrules: rename_combined_contigs

rule rename_combined_contigs:
    # standardizes header labels within contig FASTAs
    input:
        rules.combine_contigs.output.combined_contigs
    output:
        "{folder}/combined_contigs.fasta"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    params:
        prefix='C'
    shell:
        """rename.sh in={input} out={output} ow=t prefix={params.prefix}"""


rule combined_contigs_stats:
    input:
        rules.rename_combined_contigs.output
    output:
        "{folder}/combined_contigs_stats.txt"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        1
    resources:
        mem = config["mem"]
    shell:
        "stats.sh in={input} format=3 -Xmx{resources.mem}G > {output}"


rule call_genes:
    input:
        COMBINED_CONTIGS
    output:
        fna="combined/annotation/prodigal/predicted_genes.fna",
        faa="combined/annotation/prodigal/predicted_genes.faa",
        gff="combined/annotation/prodigal/predicted_genes.gff"

    conda:
        "%s/prodigal.yaml" % CONDAENV
    log:
        "logs/predict_genes.log"
    threads:
        1
    shell:
        """
            prodigal -i {input} -o {output.gff} -d {output.fna} -a {output.faa} -p meta -f gff 2> {log}
        """


rule align_reads_to_combined_contigs:
    input:
        reads=get_quality_controlled_reads,
        fasta = "{folder}/{Reference}.fasta",
    output:
        sam = temp("{folder}/sequence_alignment_{Reference}/{sample}.sam.gz"),
        unmapped= expand("{{folder}}/sequence_alignment_{{Reference}}/unmapped/{{sample}}_unmapped_{fraction}.fastq.gz",fraction= MULTIFILE_FRACTIONS)
    benchmark:
        "logs/benchmarks/sequence_alignment_{Reference}/{sample}.txt"
    params:
        input= lambda wc,input : input_params_for_bbwrap(input.reads),
        maxsites = config.get("maximum_counted_map_sites", MAXIMUM_COUNTED_MAP_SITES),
        unmapped= lambda wc,output: "outu1={0},{2} outu2={1},null".format(*output.unmapped) if PAIRED_END else "outu={0}".format(*output.unmapped),
        max_distance_between_pairs=config.get('contig_max_distance_between_pairs',CONTIG_MAX_DISTANCE_BETWEEN_PAIRS),
        paired_only = 't' if config.get("contig_map_paired_only", CONTIG_MAP_PAIRED_ONLY) else 'f',
        ambiguous = 'all' if CONTIG_COUNT_MULTI_MAPPED_READS else'best',
        min_id= config.get('contig_min_id',CONTIG_MIN_ID),
    log:
        "{folder}/logs/sequence_alignment_{Reference}/{sample}.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config["mem"]
    shell:
        """    bbwrap.sh nodisk=t \
               ref={input.fasta} \
               {params.input} \
               trimreaddescriptions=t \
               outm={output.sam} \
               {params.unmapped} \
               threads={threads} \
               pairlen={params.max_distance_between_pairs} \
               pairedonly={params.paired_only} \
               mdtag=t \
               xstag=fs \
               nmtag=t \
               sam=1.3 \
               local=t \
               ambiguous={params.ambiguous} \
               secondary=t \
               ssao=t \
               maxsites={params.maxsites} \
               -Xmx{resources.mem}G \
               2> {log}

               #max_distance_between_pairs : pairlen=32000           Set max allowed distance between paired reads.
               #(insert size)=(pairlen)+(read1 length)+(read2 length)
        """
rule pileup_combined_contigs:
    input:
        fasta = "{folder}/{Reference}.fasta",
        sam=temp("{folder}/sequence_alignment_{Reference}/{sample}.sam"),
    output:
        covstats = "{folder}/sequence_alignment_{Reference}/{sample}/{sample}_coverage.txt",
        basecov=temp("{folder}/sequence_alignment_{Reference}/{sample}/{sample}_base_coverage.txt.gz"),
        covhist= "{folder}/sequence_alignment_{Reference}/{sample}/{sample}_coverage_histogram.txt.gz",
        bincov ="{folder}/sequence_alignment_{Reference}/{sample}/{sample}_coverage_binned.txt.gz",
    params:
        pileup_secondary='t' if config.get("count_multi_mapped_reads",True) else 'f'
    benchmark:
        "logs/benchmarks/pileup_{Reference}/{sample}.txt"
    log:
        "{folder}/logs/pileup_{Reference}/{sample}.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config["mem"]
    shell:
        """
            pileup.sh ref={input.fasta} in={input.sam} threads={threads} \
            -Xmx{resources.mem}G covstats={output.covstats} \
            hist={output.covhist} basecov={output.basecov} physcov secondary={params.pileup_secondary} bincov={output.bincov} 2>> {log}
        """


rule store_bam:
    input:
        "{folder}/sequence_alignment_{Reference}/{sample}.sam.gz"
    output:
        "{folder}/sequence_alignment_{Reference}/{sample}/{sample}.bam"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config["mem"]
    shell:
        """
        reformat.sh in={input} out={output} -Xmx{resources.mem}G threads={threads}
        """







localrules: combine_coverages_of_combined_contigs, combine_bined_coverages_of_combined_contigs

rule combine_coverages_of_combined_contigs:
    input:
        covstats = expand("{folder}/sequence_alignment_{Reference}/{sample}/{sample}_coverage.txt",
            sample=SAMPLES,Reference='combined_contigs',folder=combined_contigs_folder)
    output:
        "{folder}/sequence_alignment_{Reference}/combined_median_coverage.tsv".format(Reference='combined_contigs',folder=combined_contigs_folder),
        "{folder}/sequence_alignment_{Reference}/combined_readcounts.tsv".format(Reference='combined_contigs',folder=combined_contigs_folder),
        gc_stats = "{folder}/combined_contigs_stats_gc.tsv".format(folder=combined_contigs_folder)
    run:

        import pandas as pd
        import os

        combined_cov={}
        combined_N_reads={}
        for cov_file in input:

            sample= os.path.split(cov_file)[-1].split('_')[0]
            data= pd.read_csv(cov_file,index_col=0,sep='\t')

            if cov_file == input[0]:
                data[['Length','Ref_GC']].to_csv(output.gc_stats,sep='\t')

            data.loc[data.Median_fold<0,'Median_fold']=0
            combined_cov[sample]= data.Median_fold
            combined_N_reads[sample] = data.Plus_reads+data.Minus_reads

        pd.DataFrame(combined_cov).to_csv(output[0],sep='\t')
        pd.DataFrame(combined_N_reads).to_csv(output[1],sep='\t')


rule combine_bined_coverages_of_combined_contigs:
    input:
        expand("{folder}/sequence_alignment_{Reference}/{sample}/{sample}_coverage_binned.txt.gz",
            sample=SAMPLES,Reference='combined_contigs',folder=combined_contigs_folder)
    output:
        "{folder}/sequence_alignment_{Reference}/combined_coverage_binned.tsv.gz".format(Reference='combined_contigs',folder=combined_contigs_folder),
    run:

        import pandas as pd
        import os

        binCov={}
        for cov_file in input:

            sample= os.path.split(cov_file)[-1].split('_')[0]

            binCov[sample] = pd.read_csv(cov_file,compression='gzip',sep='\t',comment='#',header=None,index_col=[0,2],usecols=[0,1,2],squeeze=True)

        binCov = pd.DataFrame(binCov)
        binCov.index.names=['Contig','Position']
        binCov.to_csv(output[0],sep='\t',compression='gzip')






#TODO parameters are not generalized
if config.get("perform_genome_binning", True):
  if config['combine_contigs_params']['binner']=='concoct':

      rule run_concoct:
          input:
              coverage= "{folder}/sequence_alignment_{Reference}/combined_median_coverage.tsv".format(Reference='combined_contigs',folder=combined_contigs_folder),
              fasta= "{folder}/{Reference}.fasta".format(Reference='combined_contigs',folder=combined_contigs_folder)
          output:
              expand("{folder}/binning/{file}",folder=combined_contigs_folder,file=['means_gt2500.csv','PCA_components_data_gt2500.csv','original_data_gt2500.csv','PCA_transformed_data_gt2500.csv','pca_means_gt2500.csv','args.txt','responsibilities.csv']),
          params:
              basename= lambda wc,output: os.path.dirname(output[0]),
              Nexpected_clusters=100,
              read_length=250,
              min_length=config.get("concoct_min_contig_length",2500),
              niterations=config.get("concoct_niterations",500)
          benchmark:
              "logs/benchmarks/binning/concoct.txt"
          log:
              "{folder}/binning/log.txt".format(folder=combined_contigs_folder)
          conda:
              "%s/concoct.yaml" % CONDAENV
          threads:
              10 # concoct uses 10 threads by default, wit for update: https://github.com/BinPro/CONCOCT/issues/177
          resources:
              mem = config["mem"]
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
  else:
      raise NotImplementedError("We don't have implemented the binning method: {}\ntry 'concoct'".format(config['combine_contigs_params']['binner']))

else:

    rule counts_genes_per_sample:
        input:
            gtf = "combined/annotation/prodigal/predicted_genes.gtf",
            bam = "{folder}/sequence_alignment_{Reference}/{sample}/{sample}.bam"
        output:
            summary = "{folder}/sequence_alignment_{Reference}/{sample}/feature_counts/{sample}_counts.txt.summary",
            counts = "{folder}/sequence_alignment_{Reference}/{sample}/feature_counts/{sample}_counts.txt"
        params:
            min_read_overlap = config.get("minimum_region_overlap", MINIMUM_REGION_OVERLAP),
            paired_only= "-B" if config.get('contig_map_paired_only',CONTIG_MAP_PAIRED_ONLY) else "",
            paired_mode = "-p" if PAIRED_END else "",
            multi_mapping = "-M --fraction" if config.get("contig_count_multi_mapped_reads",CONTIG_COUNT_MULTI_MAPPED_READS) else "--primary",
            feature_counts_allow_overlap = "-O --fraction" if config.get("feature_counts_allow_overlap", FEATURE_COUNTS_ALLOW_OVERLAP) else ""
        log:
            "{sample}/logs/counts_per_region.log"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        shell:
            """featureCounts \
                    --minOverlap {params.min_read_overlap} \
                    {params.paired_mode} \
                    {params.paired_only} \
                   -F GTF \
                   -T {threads} \
                   {params.multi_mapping} \
                   {params.feature_counts_allow_overlap} \
                   -t CDS \
                   -g ID \
                   -a {input.gtf} \
                   -o {output.counts} \
                   {input.bam} 2> {log}"""

rule combine_gene_counts:
    input:
        expand("{folder}/sequence_alignment_{Reference}/{sample}/feature_counts/{sample}_counts.txt",
            sample=SAMPLES,Reference='combined_contigs',folder=combined_contigs_folder)
    output:
        'contigs/combined_gene_counts.tsv',
        'contigs/combined_gene_info.tsv'
    run:
        import pandas as pd
        import os
        C= {}



        for file in input:
            D= pd.read_csv(file,index_col=0,comment='#',sep='\t')
            # contigs/sequence_alignment_combined_contigs/S1/S1.bam
            sample= D.columns[-1].split('/')[-2]
            C[sample]= D.iloc[:,-1]
        C= pd.DataFrame(C)
        C.to_csv(output[0],sep='\t')

        D.iloc[:,:-1].to_csv(output[1],sep='\t')





# # TODO: predict genes on all contigs
# # HACK: treat 'combined' as a sample name.
#
localrules: merge_combined_contig_tables
rule merge_combined_contig_tables:
    input:
        prokka = "{sample}/annotation/prokka/{sample}_plus.tsv".format(sample='combined'),
        refseq = "{sample}/annotation/refseq/{sample}_tax_assignments.tsv".format(sample='combined'),
        #counts = "{sample}/annotation/feature_counts/{sample}_counts.txt".format(sample='combined') # runMaxbin not suported as there is not one profile, but one per sample
    output:
        "{sample}/{sample}_annotations.txt".format(sample='combined')
    shell:
        "  atlas merge-tables \
             {input.prokka} \
             {input.refseq} \
             {output}"
