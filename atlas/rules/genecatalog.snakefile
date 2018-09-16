import os

localrules: concat_genes



rule concat_genes:
    input:
        expand("{sample}/annotation/predicted_genes/{sample}.fna", sample=SAMPLES)
    output:
        unfiltered=temp("gene_catalog/all_predicted_genes_unfiltered.fna"),
        all_genes= temp("gene_catalog/all_predicted_genes.fna"),
        lhist= "gene_catalog/stats/length_hist.txt"
    log:
        "log/gene_catalog/filter_genes.fasta"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", JAVA_MEM),
        java_mem = int(config.get("java_mem", JAVA_MEM) * JAVA_MEM_FRACTION)
    params:
        min_length=100
    shell:
        " cat {input} >  {output.unfiltered} ;"
        " reformat.sh "
        " in={output.unfiltered}"
        " minlength={params.min_length} "
        " lhist={output.lhist} "
        " ow=t out={output.all_genes} "
        " 2> {log} "




rule cluster_catalog:
    input:
        rules.concat_genes.output.all_genes
    output:
        "gene_catalog/gene_catalog.fna",
        "gene_catalog/gene_catalog.clstr"
    conda:
        "%s/cd-hit.yaml" % CONDAENV
    log:
        "logs/gene_catalog/cluster_genes.log"
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
            cd-hit-est -i {input} -T {threads} \
            -M {resources.mem}000 -o {params.prefix} \
            -c {params.identity} -n 9  -d 0 \
            -aS {params.coverage} -aL {params.coverage} &> >(tee {log})
            mv {params.prefix} {output[0]}
        """



# generalized rule so that reads from any "sample" can be aligned to contigs from "sample_contigs"
rule align_reads_to_gene_catalog:
    input:
        reads=expand("{{sample}}/sequence_quality_control/{{sample}}_QC_{fraction}.fastq.gz",
                                           fraction=MULTIFILE_FRACTIONS),
        fasta = "gene_catalog/gene_catalog.fna",
    output:
        sam = temp("gene_catalog/alignments/{sample}.sam")
    params:
        input = lambda wc, input : 'in='+','.join(input.reads),
        maxsites = 2,
        ambiguous = 'all',
        min_id = 0.95,
        maxindel = 1 # default 16000 good for genome deletions but not necessarily for alignment to contigs
    shadow:
        "shallow"
    log:
        "logs/gene_catalog/alignment/{sample}_map.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", JAVA_MEM),
        java_mem = int(config.get("java_mem", JAVA_MEM) * JAVA_MEM_FRACTION)
    shell:
        """
        bbwrap.sh nodisk=t \
            local=t \
            ref={input.fasta} \
            {params.input} \
            trimreaddescriptions=t \
            outm={output.sam} \
            threads={threads} \
            minid={params.min_id} \
            mdtag=t \
            xstag=fs \
            nmtag=t \
            sam=1.3 \
            local=t \
            ambiguous={params.ambiguous} \
            secondary=t \
            saa=f \
            maxsites={params.maxsites} \
            -Xmx{resources.java_mem}G \
            2> {log}
        """


rule pileup_gene_cluster:
    input:
        sam = "gene_catalog/alignments/{sample}.sam",
        bam = "gene_catalog/alignments/{sample}.bam"
    output:
        covstats = temp("gene_catalog/alignments/{sample}_coverage.tsv"),
    params:
        pileup_secondary = 't' # a read maay map to different genes
    log:
        "logs/gene_catalog/alignment/{sample}_pileup.log"
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
               secondary={params.pileup_secondary} \
                2> {log}
        """


localrules: combine_gene_coverages

rule combine_gene_coverages:
    input:
        covstats = expand("gene_catalog/alignments/{sample}_coverage.tsv",
            sample=SAMPLES)
    output:
        "gene_catalog/counts/raw_counts.tsv"
    run:

        import pandas as pd
        import os


        combined_N_reads={}
        for cov_file in input:

            sample= os.path.split(cov_file)[-1].split('_')[0]
            data= pd.read_table(cov_file,index_col=0)
            combined_N_reads[sample] = data.Plus_reads+data.Minus_reads

        pd.DataFrame(combined_N_reads).to_csv(output[0],sep='\t')

localrules: get_gene_catalog_annotations
rule get_gene_catalog_annotations:
    input:
        gene_catalog= 'gene_catalog/gene_catalog.fna',
        eggNOG= expand('{sample}/annotation/eggNOG.tsv',sample=SAMPLES),
        refseq= expand('{sample}/annotation/refseq/{sample}_tax_assignments.tsv',sample=SAMPLES),
        scg= expand("gene_catalog/annotation/single_copy_genes_{domain}.tsv",domain=['bacteria','archaea'])
    output:
        annotations= "gene_catalog/annotations.tsv",
    run:
        import pandas as pd

        gene_ids=[]
        with open(input.gene_catalog) as fasta_file:
            for line in fasta_file:
                if line[0]=='>':
                    gene_ids.append(line[1:].strip().split()[0])

        eggNOG=pd.DataFrame()
        for annotation_file in input.eggNOG:
            eggNOG=eggNOG.append(pd.read_table(annotation_file, index_col=0))

        refseq=pd.DataFrame()
        for annotation_file in input.refseq:
            refseq=refseq.append(pd.read_table(annotation_file, index_col=1))

        scg=pd.DataFrame()
        for annotation_file in input.scg:
            d= pd.read_table(annotation_file, index_col=0,header=None)
            d.columns = 'scg_'+ os.path.splitext(annotation_file)[0].split('_')[-1] # bacteria or archaea
            scg=scg.append(d)


        annotations= refseq.join(eggNOG).join(scg).loc[gene_ids]
        annotations.to_csv(output.annotations,sep='\t')


rule predict_single_copy_genes:
    input:
        "gene_catalog/gene_catalog.fna"
    output:
        "gene_catalog/annotation/single_copy_genes_{domain}.tsv",
    params:
        script_dir = os.path.dirname(os.path.abspath(workflow.snakefile)),
        key = lambda wc: wc.domain[:3] #bac for bacteria, #archaea
    conda:
        "%s/DASTool.yaml" % CONDAENV # needs pearl
    threads:
        config['threads']
    run:
        " DIR=$(dirname $(which DAS_Tool)) "
        ";"
        " {params.script_dir}/rules/scg_blank_diamond.rb diamond"
        " {input} "
        " $DIR\/db/{params.key}.all.faa "
        " $DIR\/db/{params.key}.scg.faa "
        " $DIR\/db/{params.key}.scg.lookup "
        " {threads} "
        " &> >(tee {log}) "
        " mv {input[0]}.{wildcards.domain}.scg > {output}"




#
# ############## Canopy clustering
#
# rule reformat_for_canopy:
#         input:
#             "mapresults/Gene_catalog_CE/combined_Nmaped_reads.tsv"
#         output:
#             "mapresults/Gene_catalog_CE/nseq.tsv"
#         run:
#             import pandas as pd
#
#             D= pd.read_table(input[0], index_col=0)
#             D.index= D.index.map(lambda s: s.split()[0])
#             D=D.astype(int)
#             D.to_csv(output[0],sep='\t',header=False)
#
#
# rule canopy_clustering:
#     input:
#         rules.reformat_for_canopy.output
#     output:
#         cluster="mapresults/Gene_catalog_CE/canopy_cluster.tsv",
#         profile="mapresults/Gene_catalog_CE/cluster_profiles.tsv",
#     params:
#         canopy_params=config.get("canopy_params","")
#     log:
#         "mapresults/Gene_catalog_CE/canopy.log"
#     benchmark:
#         "logs/benchmarks/canopy_clustering.txt"
#     conda:
#         "%s/canopy.yaml" % CONDAENV
#     threads:
#         12
#     resources:
#         mem= 220
#     shell:
#         """
#         canopy -i {input} -o {output.cluster} -c {output.profile} -n {threads} --canopy_size_stats_file {log} {params.canopy_params} 2> >(tee {log})
#
#         """
