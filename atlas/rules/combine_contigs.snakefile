import os
import re
import sys
import tempfile
import warnings

from default_values import *

combined_contigs_folder='combined'
COMBINED_CONTIGS= "{folder}/combined_contigs.fasta".format(folder=combined_contigs_folder)

rule combine_contigs_report:
    input:
        combined_contigs= COMBINED_CONTIGS,
        combined_contigs_stats="{folder}/combined_contigs_stats.txt".format(folder=combined_contigs_folder),
        annotation= "{sample}/{sample}_annotations.txt".format(sample='combined'),
        median_coverage="{folder}/sequence_alignment_{Reference}/combined_median_coverage.tsv".format(Reference='combined_contigs',folder=combined_contigs_folder),
        gc_stats = "{folder}/combined_contigs_stats_gc.tsv".format(folder=combined_contigs_folder),
        binned_coverage = "{folder}/sequence_alignment_{Reference}/combined_coverage_binned.tsv.gz".format(Reference='combined_contigs',folder=combined_contigs_folder),
        #concoct="{folder}/binning/{file}".format(folder=combined_contigs_folder,file='means_gt2500.csv')
        #bam=expand("{folder}/sequence_alignment_{Reference}/{sample}/{sample}.bam",sample=SAMPLES,Reference='combined_contigs',folder=combined_contigs_folder),
    output:
        touch("Combined_contigs_done")

#### combine contigs

#TODO: put in conf values
config['combine_contigs']=dict(min_overlap = 200,
                               max_substitutions=4,
                               dont_allow_N=True,
                               remove_cycles=True,
                               trim_contradictions=True, #False
                               binner='concoct'
                               )

config['perform_genome_binning']= False

rule combine_contigs:
    input:
        expand("{sample}/{sample}_contigs.fa",sample=SAMPLES)
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
        mem = config.get("java_mem", JAVA_MEM)
    params:
       input=lambda wc,input: ','.join(input),
       min_length=config.get("minimum_contig_length", MINIMUM_CONTIG_LENGTH),
       min_overlap=config['combine_contigs']['min_overlap'],
       max_substitutions=config['combine_contigs']['max_substitutions'],
       dont_allow_N= 't' if config['combine_contigs']['dont_allow_N'] else 'f',
       remove_cycles='t' if config['combine_contigs']['remove_cycles'] else 'f',
       trim_contradictions='t' if config['combine_contigs']['trim_contradictions'] else 'f'

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
            -Xmx{resources.mem}G 2> >(tee {log})
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
        mem = config.get("java_mem", JAVA_MEM)
    shell:
        "stats.sh in={input} format=3 -Xmx{resources.mem}G > {output}"


# rule align_reads_to_combined_contigs:
#     input:
#         unpack(get_quality_controlled_reads),
#         fasta = "{folder}/{Reference}.fasta",
#     output:
#         sam = temp("{folder}/sequence_alignment_{Reference}/{sample}.sam"),
#         unmapped= expand("{{folder}}/sequence_alignment_{{Reference}}/unmapped/{{sample}}_unmapped_{fraction}.fastq.gz",fraction=multifile_fractions)
#     benchmark:
#         "logs/benchmarks/sequence_alignment_{Reference}/{sample}.txt"
#     params:
#         input= lambda wc,input : input_params_for_bbwrap(wc,input),
#         maxsites = config.get("maximum_counted_map_sites", MAXIMUM_COUNTED_MAP_SITES),
#         unmapped= lambda wc,output: "outu1={0},{2} outu2={1},null".format(*output.unmapped) if paired_end else "outu={0}".format(*output.unmapped),
#         max_distance_between_pairs=config.get('contig_max_distance_between_pairs',CONTIG_MAX_DISTANCE_BETWEEN_PAIRS),
#         paired_only = 't' if config.get("contig_map_paired_only", CONTIG_MAP_PAIRED_ONLY) else 'f',
#         ambiguous = 'all' if CONTIG_COUNT_MULTI_MAPPED_READS else'best',
#         min_id= config.get('contig_min_id',CONTIG_MIN_ID),
#     log:
#         "{folder}/logs/sequence_alignment_{Reference}/{sample}.log"
#     conda:
#         "%s/required_packages.yaml" % CONDAENV
#     threads:
#         config.get("threads", 1)
#     resources:
#         mem = config.get("java_mem", JAVA_MEM)
#     shell:
#         """{SHPFXM} bbwrap.sh nodisk=t \
#                ref={input.fasta} \
#                {params.input} \
#                trimreaddescriptions=t \
#                outm={output.sam} \
#                {params.unmapped} \
#                threads={threads} \
#                pairlen={params.max_distance_between_pairs} \
#                pairedonly={params.paired_only} \
#                mdtag=t \
#                xstag=fs \
#                nmtag=t \
#                sam=1.3 \
#                local=t \
#                ambiguous={params.ambiguous} \
#                secondary=t \
#                ssao=t \
#                maxsites={params.maxsites} \
#                -Xmx{resources.mem}G \
#                2> {log}
#
#                #max_distance_between_pairs : pairlen=32000           Set max allowed distance between paired reads.
#                #(insert size)=(pairlen)+(read1 length)+(read2 length)
#         """
# rule pileup_combined_contigs:
#     input:
#         fasta = "{folder}/{Reference}.fasta",
#         sam=temp("{folder}/sequence_alignment_{Reference}/{sample}.sam"),
#     output:
#         covstats = "{folder}/sequence_alignment_{Reference}/{sample}/{sample}_coverage.txt",
#         basecov=temp("{folder}/sequence_alignment_{Reference}/{sample}/{sample}_base_coverage.txt.gz"),
#         covhist= "{folder}/sequence_alignment_{Reference}/{sample}/{sample}_coverage_histogram.txt.gz",
#         bincov ="{folder}/sequence_alignment_{Reference}/{sample}/{sample}_coverage_binned.txt.gz",
#     params:
#         pileup_secondary='t' if config.get("count_multi_mapped_reads",True) else 'f'
#     benchmark:
#         "logs/benchmarks/pileup_{Reference}/{sample}.txt"
#     log:
#         "{folder}/logs/pileup_{Reference}/{sample}.log"
#     conda:
#         "%s/required_packages.yaml" % CONDAENV
#     threads:
#         config.get("threads", 1)
#     resources:
#         mem = config.get("java_mem", JAVA_MEM)
#     shell:
#         """
#             pileup.sh ref={input.fasta} in={input.sam} threads={threads} \
#             -Xmx{resources.mem}G covstats={output.covstats} \
#             hist={output.covhist} basecov={output.basecov} physcov secondary={params.pileup_secondary} bincov={output.bincov} 2>> {log}
#         """


# rule store_bam:
#     input:
#         "{folder}/sequence_alignment_{Reference}/{sample}.sam.gz"
#     output:
#         "{folder}/sequence_alignment_{Reference}/{sample}/{sample}.bam"
#     conda:
#         "%s/required_packages.yaml" % CONDAENV
#     threads:
#         config.get("threads", 1)
#     resources:
#         mem = config.get("java_mem", JAVA_MEM)
#     shell:
#         """
#         reformat.sh in={input} out={output} -Xmx{resources.mem}G threads={threads}
#         """







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
            data= pd.read_table(cov_file,index_col=0)

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

            binCov[sample] = pd.read_table(cov_file,compression='gzip',comment='#',header=None,index_col=[0,2],usecols=[0,1,2],squeeze=True)

        binCov = pd.DataFrame(binCov)
        binCov.index.names=['Contig','Position']
        binCov.to_csv(output[0],sep='\t',compression='gzip')



#TODO parameters are not generalized
if config.get("perform_genome_binning", True):
  if config['combine_contigs']['binner']=='concoct':

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
  else:
      raise NotImplementedError("We don't have implemented the binning method: {}\ntry 'concoct'".format(config['combine_contigs']['binner']))

else:
    rule run_prokka_annotation:
        input:
            "{folder}/sequence_alignment_{Reference}/combined_median_coverage.tsv".format(Reference='combined_contigs',folder=combined_contigs_folder)
        output:
            discrepancy = "{sample}/annotation/prokka/{sample}.err",
            faa = "{sample}/annotation/prokka/{sample}.faa",
            ffn = "{sample}/annotation/prokka/{sample}.ffn",
            fna = "{sample}/annotation/prokka/{sample}.fna",
            fsa = "{sample}/annotation/prokka/{sample}.fsa",
            gbk = "{sample}/annotation/prokka/{sample}.gbk",
            gff = "{sample}/annotation/prokka/{sample}.gff",
            log = "{sample}/annotation/prokka/{sample}.log",
            sqn = "{sample}/annotation/prokka/{sample}.sqn",
            tbl = "{sample}/annotation/prokka/{sample}.tbl",
            tsv = "{sample}/annotation/prokka/{sample}.tsv",
            txt = "{sample}/annotation/prokka/{sample}.txt"
        benchmark:
            "logs/benchmarks/prokka/{sample}.txt"
        params:
            outdir = lambda wc, output: os.path.dirname(output.faa),
            kingdom = config.get("prokka_kingdom", PROKKA_KINGDOM)
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        shell:
            """prokka --outdir {params.outdir} \
                   --force \
                   --prefix {wildcards.sample} \
                   --locustag {wildcards.sample} \
                   --kingdom {params.kingdom} \
                   --metagenome \
                   --cpus {threads} \
                   {input}"""


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
