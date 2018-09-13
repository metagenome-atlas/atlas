import os

localrules: concat_genes,rename_combined_contigs


rule combine_contigs:
    input:
        expand("{sample}/assembly/{sample}_prefilter_contigs.fasta",sample=SAMPLES)
    output:
        combined_contigs=temp("{folder}/combined_contigs_oldnames.fasta"),
        cluster_stats="{folder}/combined_contigs_kmerfrequencies.txt",
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
       min_overlap=200,
       max_substitutions=4,
       dont_allow_N='t',
       remove_cycles='t',
       trim_contradictions='f'



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

#dot -Tpdf combined_contigs_graph.dot -o combined_clusters.pdf

rule rename_combined_contigs:
    # standardizes header labels within contig FASTAs
    input:
        "{folder}/combined_contigs_oldnames.fasta"
    output:
        "{folder}/combined_contigs.fasta"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    params:
        prefix='C'
    shell:
        """rename.sh in={input} out={output} ow=t prefix={params.prefix}"""


rule gene_calling_combined_contigs:
    input:
        "{folder}/{Reference}.fasta"
    output:
        fna="{folder}/predicted_genes/{Reference}_genes.fna",
        faa="{folder}/predicted_genes/{Reference}_genes.faa",
        gff="{folder}/predicted_genes/{Reference}_genes.gff"

    conda:
        "%s/gene_catalog.yaml" % CONDAENV
    log:
        "{folder}/logs/predict_genes_{Reference}.log"
    threads:
        1
    shell:
        """
            prodigal -i {input} -o {output.gff} -d {output.fna} -a {output.faa} -p meta -f gff 2> >(tee {log})
        """

# TODO:  max insertsize per sample
# Max insertzize and readlen per sample and also for combined. !!
rule align_reads_to_combined_contigs:
    input:
        unpack(get_quality_controlled_reads),
        fasta = "{folder}/{Reference}.fasta",
    output:
        sam = temp("{folder}/sequence_alignment_{Reference}/{sample}.sam.gz"),
        unmapped= expand("{{folder}}/sequence_alignment_{{Reference}}/unmapped/{{sample}}_unmapped_{fraction}.fastq.gz",fraction=multifile_fractions)
    benchmark:
        "logs/benchmarks/sequence_alignment_{Reference}/{sample}.txt"
    params:
        input= lambda wc,input : input_params_for_bbwrap(wc,input),
        maxsites = config.get("maximum_counted_map_sites", MAXIMUM_COUNTED_MAP_SITES),
        unmapped= lambda wc,output: "outu1={0},{2} outu2={1},null".format(*output.unmapped) if paired_end else "outu={0}".format(*output.unmapped),
        max_distance_between_pairs=1000,
        samestrandpairs='f',
        paired_only='t',
        ambiguous='best'
    log:
        "{folder}/logs/sequence_alignment_{Reference}/{sample}.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", JAVA_MEM)
    shell:
        """{SHPFXM} bbwrap.sh nodisk=t \
               ref={input.fasta} \
               {params.input} \
               trimreaddescriptions=t \
               outm={output.sam} \
               {params.unmapped} \
               threads={threads} \
               pairlen={params.max_distance_between_pairs} \
               pairedonly={params.paired_only} \
               samestrandpairs={params.samestrandpairs} \
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
rule pileup:
    input:
        fasta = "{folder}/{Reference}.fasta",
        sam="{folder}/sequence_alignment_{Reference}/{sample}.sam.gz",
        predicted_genes="{folder}/predicted_genes/{Reference}_genes.fna"
    output:
        covstats = "{folder}/sequence_alignment_{Reference}/{sample}/{sample}_coverage.txt",
        gene_coverage="{folder}/sequence_alignment_{Reference}/{sample}/{sample}_gene_coverage.txt.gz",
        basecov="{folder}/sequence_alignment_{Reference}/{sample}/{sample}_base_coverage.txt.gz",
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
        mem = config.get("java_mem", JAVA_MEM)
    shell:
        """
            pileup.sh ref={input.fasta} in={input.sam} threads={threads} \
            -Xmx{resources.mem}G covstats={output.covstats} \
            fastaorf={input.predicted_genes} outorf={output.gene_coverage} \
            hist={output.covhist} basecov={output.basecov} physcov secondary={params.pileup_secondary} bincov={output.bincov} 2>> {log}
        """

localrules: combine_coverages_of_combined_contigs

rule combine_coverages_of_combined_contigs:
    input:
        covstats = expand("{folder}/sequence_alignment_{Reference}/{sample}/{sample}_coverage.txt",
            sample=SAMPLES,Reference='combined_contigs',folder='combined')
    output:
        "{folder}/sequence_alignment_{Reference}/combined_median_coverage.tsv".format(Reference='combined_contigs',folder='combined'),
        "{folder}/sequence_alignment_{Reference}/combined_readcounts.tsv".format(Reference='combined_contigs',folder='combined')
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

localrules: run_concoct
#TODO parameters are not generalized

rule run_concoct:
    input:
        coverage= "{folder}/sequence_alignment_{Reference}/combined_median_coverage.tsv".format(Reference='combined_contigs',folder='combined'),
        fasta= "{folder}/{Reference}.fasta".format(Reference='combined_contigs',folder='combined')
    output:
        "{folder}/binning/combined_readcounts.tsv".format(folder='combined'),
    params:
        basename= lambda wc,output: os.path.dirname(output[0]),
        Nexpected_clusters=100,
        read_length=250,
        min_length=config.get("concoct_min_contig_length",2500),
        niterations=config.get("concoct_niterations",500)
    benchmark:
        "logs/benchmarks/binning/concoct.txt"
    log:
        "{folder}/binning/log.txt".format(folder='combined')
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



# alredy possible with prefilter contig stats
# rule calculate_combined_contigs_stats:
#     input:
#         rules.rename_combined_contigs.output
#     output:
#         "{folder}/{Reference}_stats.txt"
#     conda:
#         "%s/required_packages.yaml" % CONDAENV
#     threads:
#         1
#     resources:
#         mem = config.get("java_mem", JAVA_MEM)
#     shell:
#         "{SHPFXS} stats.sh in={input} format=3 -Xmx{resources.mem}G > {output}"



# rule copy_contigs_for_annotation:
#     input:
#         expand(rules.rename_combined_contigs.output,folder='combined')
#     output:
#         "{sample}/{sample}_contigs.fasta".format(sample='combined')
#     shell:
#         """
#             cp {input} {output}
#         """

rule merge_combined_contig_tables:
    input:
        prokka = "{sample}/annotation/prokka/{sample}_plus.tsv".format(sample='combined'),
        refseq = "{sample}/annotation/refseq/{sample}_tax_assignments.tsv".format(sample='combined'),
        #counts = "{sample}/annotation/feature_counts/{sample}_counts.txt"
    output:
        "{sample}/{sample}_annotations.txt".format(sample='combined')
    shell:
        "  atlas merge-tables \
             {input.prokka} \
             {input.refseq} \
             {output}"


####### make gene catalog


rule gene_calling:
    input:
        "{sample}/assembly/{sample}_prefilter_contigs.fasta"
    output:
        fna="{sample}/assembly/predicted_genes.fna",
        faa="{sample}/assembly/predicted_genes.faa",
        gff="{sample}/assembly/predicted_genes.gff"

    conda:
        "%s/gene_catalog.yaml" % CONDAENV
    log:
        "{sample}/logs/{sample}_gene_calling.log"
    threads:
        1
    shell:
        """
            prodigal -i {input} -o {output.gff} -d {output.fna} -a {output.faa} -p meta -f gff 2> >(tee {log})
        """

rule concat_genes:
    input:
        expand("{sample}/assembly/predicted_genes.faa", sample=SAMPLES)
    output:
        temp("gene_catalog/all_predicted_genes.faa")
    shell:
        """
            cat {input} >  {output}
        """

rule cluster_catalog:
    input:
        rules.concat_genes.output
    output:
        "gene_catalog/gene_catalog.fasta",
        "gene_catalog/gene_catalog.clstr"
    conda:
        "%s/gene_catalog.yaml" % CONDAENV
    log:
        "logs/make_unique_gene_catalog.log"
    threads:
        8
    resources:
        mem=20
    params:
        prefix= lambda wc,output: os.path.splitext(output[0])[0],
        coverage=0.9,
        identity=0.95
    shell:
        """
            cd-hit-est -i {input} -T {threads} -M {resources.mem}000 -o {params.prefix} -c {params.identity} -n 9  -d 0 -aS {params.coverage} -aL {params.coverage} 2> >(tee {log})
            mv {params.prefix} {output[0]}
        """

############## Canopy clustering

rule reformat_for_canopy:
        input:
            "mapresults/Gene_catalog_CE/combined_Nmaped_reads.tsv"
        output:
            "mapresults/Gene_catalog_CE/nseq.tsv"
        run:
            import pandas as pd

            D= pd.read_table(input[0], index_col=0)
            D.index= D.index.map(lambda s: s.split()[0])
            D=D.astype(int)
            D.to_csv(output[0],sep='\t',header=False)


rule canopy_clustering:
    input:
        rules.reformat_for_canopy.output
    output:
        cluster="mapresults/Gene_catalog_CE/canopy_cluster.tsv",
        profile="mapresults/Gene_catalog_CE/cluster_profiles.tsv",
    params:
        canopy_params=config.get("canopy_params","")
    log:
        "mapresults/Gene_catalog_CE/canopy.log"
    benchmark:
        "logs/benchmarks/canopy_clustering.txt"
    conda:
        "%s/canopy.yaml" % CONDAENV
    threads:
        12
    resources:
        mem= 220
    shell:
        """
        canopy -i {input} -o {output.cluster} -c {output.profile} -n {threads} --canopy_size_stats_file {log} {params.canopy_params} 2> >(tee {log})

        """
