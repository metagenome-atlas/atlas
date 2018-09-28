import os


rule gene_catalog:
    input:
        "Genecatalog/protein_catalog.faa",
        "Genecatalog/gene_catalog.fna",
        "Genecatalog/gene_catalog.faa",
        "Genecatalog/counts/median_coverage.tsv",
        expand("Genecatalog/annotation/single_copy_genes_{domain}.tsv",domain=['bacteria','archaea'])


localrules: concat_genes
rule concat_genes:
    input:
        faa= expand("{sample}/annotation/predicted_genes/{sample}.faa", sample=SAMPLES),
        fna= expand("{sample}/annotation/predicted_genes/{sample}.fna", sample=SAMPLES)
    output:
        faa=  temp("Genecatalog/all_genes_unfiltered.faa"),
        fna = temp("Genecatalog/all_genes_unfiltered.fna"),
    shell:
        " cat {input.faa} >  {output.faa} ;"
        " cat {input.fna} > {output.fna}"



localrules: filter_genes
rule filter_genes:
    input:
        fna="Genecatalog/all_genes_unfiltered.fna",
        faa="Genecatalog/all_genes_unfiltered.fna"
    output:
        fna= "Genecatalog/all_genes/predicted_genes.fna",
        faa= "Genecatalog/all_genes/predicted_genes.faa",
    threads:
        1
    params:
        min_length=config['genecatalog']['minlength_nt']
    run:
        from Bio import SeqIO
        faa = SeqIO.parse(input.faa,'fasta')
        fna = SeqIO.parse(input.fna,'fasta')

        with open(output.faa,'w') as out_faa, open(output.fna,'w') as out_fna:

            for gene in fna:
                protein = next(faa)

                if len(gene) >= params.min_length:
                    SeqIO.write(gene,out_fna,'fasta')
                    SeqIO.write(protein,out_faa,'fasta')


rule cluster_proteins:
    input:
        faa= "Genecatalog/all_genes/predicted_genes.faa"
    output:
        db=temp(expand("Genecatalog/all_genes/predicted_genes.{dbext}",dbext=['db','db.dbtype',
                                                                        'db.index', 'db.lookup', 'db_h', 'db_h.index'])),
        clusterdb = temp(expand("Genecatalog/clustering/protein_clusters.{ext}",ext=['db','db.index'])),
    conda:
        "%s/mmseqs.yaml" % CONDAENV
    log:
        "logs/Genecatalog/clustering/cluster_proteins.log"
    threads:
        config.get("threads", 1)
    params:
        tmpdir= temp(directory(os.path.join(config['tmpdir'],"mmseqs"))),
        clustermethod = config['cluster_proteins']['method'],
        coverage=config['cluster_proteins']['coverage'], #0.8,
        evalue=config['cluster_proteins']['evalue'], # 0.001
        minid=config['cluster_proteins']['minid'], # 0.00
        extra=config['cluster_proteins']['extra']
    shell:
        """
            mmseqs createdb {input.faa} {output.db[0]} > >(tee  {log})

            mkdir -p {params.tmpdir}

            mmseqs {params.clustermethod} -c {params.coverage} -e {params.evalue} \
            --min-seq-id {params.minid} {params.extra} \
            --threads {threads} {output.db[0]} {output.clusterdb[0]} {params.tmpdir}  > >(tee -a  {log})
        """


rule get_rep_proteins:
    input:
        db= rules.cluster_proteins.output.db,
        clusterdb = rules.cluster_proteins.output.clusterdb,
    output:
        cluster_attribution = temp("Genecatalog/orf2proteins_oldnames.tsv"),
        rep_seqs_db = temp(expand("Genecatalog/protein_catalog.{exp}",exp=['db','db.index'])),
        rep_seqs = temp("Genecatalog/protein_catalog_oldnames.faa")
    conda:
        "%s/mmseqs.yaml" % CONDAENV
    log:
        "logs/Genecatalog/clustering/get_rep_proteins.log"
    threads:
        config.get("threads", 1)
    shell:
        """
        mmseqs createtsv {input.db[0]} {input.db[0]} {input.clusterdb[0]} {output.cluster_attribution}  > >(tee   {log})

        mmseqs result2repseq {input.db[0]} {input.clusterdb[0]} {output.rep_seqs_db[0]}  > >(tee -a  {log})
        mmseqs result2flat {input.db[0]} {input.db[0]} {output.rep_seqs_db[0]} {output.rep_seqs}  > >(tee -a  {log})

        """

def gen_names_for_range(N,prefix='',start=1):
    """generates a range of IDS with leading zeros so sorting will be ok"""
    n_leading_zeros= len(str(N))
    format_int=prefix+'{:0'+str(n_leading_zeros)+'d}'
    return [format_int.format(i) for i in range(start,N+start)]



localrules: rename_protein_catalog
rule rename_protein_catalog:
    input:
        cluster_attribution = "Genecatalog/orf2proteins_oldnames.tsv",
        rep_seqs = "Genecatalog/protein_catalog_oldnames.faa"
    output:
        cluster_attribution = "Genecatalog/clustering/orf2proteins.tsv",
        rep_seqs = "Genecatalog/protein_catalog.faa"
    run:
        import pandas as pd
        # CLuterID    GeneID    empty third column
        gene2proteins= pd.read_table(input.cluster_attribution,index_col=1, header=None)

        protein_clusters_old_names= gene2proteins[0].unique()

        map_names = dict(zip(protein_clusters_old_names,
                             gen_names_for_range(len(protein_clusters_old_names),'ProteinCluster')))

        gene2proteins['proteinCluster'] = gene2proteins[0].map(map_names)
        gene2proteins.index.name='ORF'
        gene2proteins['proteinCluster'].to_csv(output.cluster_attribution,sep='\t',header=True)


        with open(output.rep_seqs,'w') as fout:
            with open(input.rep_seqs) as fin :
                for line in fin:
                    if line[0]=='>':
                        fout.write(">{new_name}\n".format(new_name=map_names[line[1:].strip()]))
                    else:
                        fout.write(line)

# cluster genes


# Create fasta files of the genes (nt) for all genes in a protein cluster for subclustering
rule dispatch_fasta:
    input:
        genes2proteins = "Genecatalog/clustering/orf2proteins.tsv",
        fna= "Genecatalog/all_genes/predicted_genes.fna",
    output:
        temp_folder= temp(directory("Genecatalog/clustering/nt_unclustered")),
        unique_fna = temp("Genecatalog/clustering/unique_genes.fna")
    run:
        import pandas as pd
        import os, shutil
        from Bio import SeqIO

        os.mkdir(output.temp_folder)

        genes2proteins=pd.read_table(input.genes2proteins,
                             index_col=0,squeeze=True)
        #no need to cluster unique genes
        unique_genes= genes2proteins.drop_duplicates(keep=False)

        if unique_genes.shape[0]>0:
            genes2proteins.loc[unique_genes.index]='unique'
        else:
            open(output.unique_fna,'w').close()

        # create individual file handles
        from collections import defaultdict
        output_groups = defaultdict([])
        #{dict( (proteinID,open(os.path.join(output.temp_folder, f"{proteinID}.fna"),"w"))
                          #for proteinID in genes2proteins.unique())}

        for seq in SeqIO.parse(input.fna,'fasta'):
            output_groups[genes2proteins[seq.name]].append(seq)

        for proteinID in output_groups:
            SeqIO.write(output_groups[proteinID],f"{output.temp_folder}/{proteinID}.fna",'fasta')

        shutil.move(os.path.join(output.temp_folder, "unique.fna"),output.unique_fna)


rule subcluster_genes:
    input:
        fna_dir=directory("Genecatalog/clustering/nt_unclustered"),
    output:
        temp(directory("Genecatalog/clustering/gene_subclusters"))
    conda:
        "%s/cd-hit.yaml" % CONDAENV
    log:
        "logs/Genecatalog/cluster_genes.log"
    threads:
        config.get("threads", 1)
    resources:
        mem=config.get("java_mem", JAVA_MEM)
    params:
        snakefile= os.path.join(os.path.dirname(workflow.snakefile),'rules','cluster_genes.snakefile'),
        coverage=config['genecatalog']['coverage'],
        identity=config['genecatalog']['minid']
    shell:
        " snakemake -s {params.snakefile} "
        " --config  unclustered_dir={input.fna_dir} clustered_dir={output} "
        " identity={params.identity} coverage={params.coverage} "
        " -j {threads} -p "
        " --resources mem={resources.mem} "
        " &> {log} "





localrules: combine_gene_clusters, parse_clstr_files, rename_gene_clusters
rule combine_gene_clusters:
    input:
        clustered_dir= directory("Genecatalog/clustering/gene_subclusters"),
        unique_fna = "Genecatalog/clustering/unique_genes.fna",
    output:
        temp("Genecatalog/gene_catalog_oldnames.fna"),
    shell:
        "cat {input.clustered_dir}/*.fna {input.unique_fna} > {output}"

def parse_cd_hit_file(clstr_file):
    """

    >Cluster 0
    0	342nt, >S1_83_1... *
    1	342nt, >S2_82_1... at +/100.00%
    >Cluster 1
    0	339nt, >S1_61_1... *
    1	339nt, >S2_59_1... at +/100.00%


    """

    def parse_line(line):
        _, length, name, identity = line.strip().replace('...','\t').replace(', ','\t').split('\t')

        length= int(length.replace('nt',''))
        name=name[1:]
        if '*' in identity:
            identity= None
        else:
            identity= float(identity[identity.rfind('/')+1:identity.rfind('%')])

        return name,length, identity

    Clusters= []
    with open(clstr_file) as f:
        for line in f:
            if line[0]=='>': #new cluster
                cluster= dict(elements=[],representative=None)
                Clusters.append(cluster)
            else:
                name,length, identity =parse_line(line)
                cluster['elements'].append((name,length, identity))
                if identity is None:
                    cluster['representative']= name
    return Clusters


def write_cd_hit_clusters(Clusters,file_handle):
        for cluster in Clusters:
            for element in cluster['elements']:
                file_handle.write(f"{element[0]}\t{element[1]}\t{element[2]}\t{cluster['representative']}\n")



rule parse_clstr_files:
    input:
        clustered_dir= directory("Genecatalog/clustering/gene_subclusters")
    output:
        temp("Genecatalog/clustering/orf2gene_oldnames.tsv")
    run:
        from glob import glob

        with open(output[0],'w') as fout:
            fout.write(f"ORF\tLength\tIdentity\tGene\n")
            for clstr_file in glob("{input.clustered_dir}/*.clstr"):
                Clusters=parse_cd_hit_file(clstr_file)
                write_cd_hit_clusters(Clusters,fout)



rule rename_gene_clusters:
    input:
        fna = "Genecatalog/gene_catalog_oldnames.fna",
        faa= "Genecatalog/all_genes/predicted_genes.faa",
        orf2gene = "Genecatalog/clustering/orf2gene_oldnames.tsv",
        orf2protein = "Genecatalog/orf2proteins_oldnames.tsv"
    output:
        fna= "Genecatalog/gene_catalog.fna",
        faa= "Genecatalog/gene_catalog.faa",
        mapping_file = "Genecatalog/clustering/gene_mapping_file.tsv"
    run:
        import pandas as pd
        from Bio import SeqIO

        orf2gene= pd.read_table(input.orf2gene,index_col=0)
        orf2protein = pd.read_table(input.orf2protein,index_col=0)

        orf2gene=orf2gene.join(orf2gene)

        # rename gene repr to Gene0000XX

        gene_clusters_old_names= gene2proteins['Gene'].unique()

        map_names = dict(zip(gene_clusters_old_names,
                             gen_names_for_range(len(gene_clusters_old_names),'Gene')))

        gene2proteins['Gene'] = gene2proteins['Gene'](map_names)
        gene2proteins.to_csv(output.mapping_file,sep='\t',header=True)

        # rename fna
        with open(output.fna,'w') as f_out:
            for gene in SeqIO.parse(output.fna,'fasta'):
                gene.id = old2new_names[gene.name]
                SeqIO.write(gene,f_out,'fasta')

        # rename faa
        with open(output.faa,'w') as f_out:
            for gene in SeqIO.parse(input.faa,'fasta'):
                if gene.name in old2new_names:
                    gene.id = old2new_names[gene.name]
                    SeqIO.write(gene,f_out,'fasta')





rule align_reads_to_Genecatalog:
    input:
        unpack(get_quality_controlled_reads),
        fasta = "Genecatalog/gene_catalog.fna",
    output:
        sam = temp("Genecatalog/alignments/{sample}.sam")
    params:
        input = lambda wc, input : input_params_for_bbwrap(wc, input),
        maxsites = 2,
        ambiguous = 'all',
        minid = config['genecatalog']['minid'],
        maxindel = 1 # default 16000 good for genome deletions but not necessarily for alignment to contigs
    shadow:
        "shallow"
    log:
        "logs/Genecatalog/alignment/{sample}_map.log"
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
            minid={params.minid} \
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


rule pileup_Genecatalog:
    input:
        sam = "Genecatalog/alignments/{sample}.sam",
        bam = "Genecatalog/alignments/{sample}.bam"
    output:
        covstats = temp("Genecatalog/alignments/{sample}_coverage.tsv"),
        basecov = temp("Genecatalog/alignments/{sample}_base_coverage.txt.gz"),
    params:
        pileup_secondary = 't' # a read maay map to different genes
    log:
        "logs/Genecatalog/alignment/{sample}_pileup.log"
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
               basecov={output.basecov} \
               secondary={params.pileup_secondary} \
                2> {log}
        """

localrules: combine_gene_coverages
rule combine_gene_coverages:
    input:
        covstats = expand("Genecatalog/alignments/{sample}_coverage.tsv",
            sample=SAMPLES)
    output:
        "Genecatalog/counts/median_coverage.tsv",
        "Genecatalog/counts/Nmapped_reads.tsv",
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


#
# localrules: get_Genecatalog_annotations
# rule get_Genecatalog_annotations:
#     input:
#         Genecatalog= 'Genecatalog/gene_catalog.fna".fna',
#         eggNOG= expand('{sample}/annotation/eggNOG.tsv',sample=SAMPLES),
#         refseq= expand('{sample}/annotation/refseq/{sample}_tax_assignments.tsv',sample=SAMPLES),
#         scg= expand("Genecatalog/annotation/single_copy_genes_{domain}.tsv",domain=['bacteria','archaea'])
#     output:
#         annotations= "Genecatalog/annotations.tsv",
#     run:
#         import pandas as pd
#
#         gene_ids=[]
#         with open(input.Genecatalog) as fasta_file:
#             for line in fasta_file:
#                 if line[0]=='>':
#                     gene_ids.append(line[1:].strip().split()[0])
#
#         eggNOG=pd.DataFrame()
#         for annotation_file in input.eggNOG:
#             eggNOG=eggNOG.append(pd.read_table(annotation_file, index_col=0))
#
#         refseq=pd.DataFrame()
#         for annotation_file in input.refseq:
#             refseq=refseq.append(pd.read_table(annotation_file, index_col=1))
#
#         scg=pd.DataFrame()
#         for annotation_file in input.scg:
#             d= pd.read_table(annotation_file, index_col=0,header=None)
#             d.columns = 'scg_'+ os.path.splitext(annotation_file)[0].split('_')[-1] # bacteria or archaea
#             scg=scg.append(d)
#
#
#         annotations= refseq.join(eggNOG).join(scg).loc[gene_ids]
#         annotations.to_csv(output.annotations,sep='\t')



rule predict_single_copy_genes:
    input:
        "Genecatalog/gene_catalog.fna"
    output:
        "Genecatalog/annotation/single_copy_genes_{domain}.tsv",
    params:
        script_dir = os.path.dirname(os.path.abspath(workflow.snakefile)),
        key = lambda wc: wc.domain[:3] #bac for bacteria, #archaea
    conda:
        "%s/DASTool.yaml" % CONDAENV # needs pearl
    threads:
        config['threads']
    shell:
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
#             "mapresults/Genecatalog_CE/combined_Nmaped_reads.tsv"
#         output:
#             "mapresults/Genecatalog_CE/nseq.tsv"
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
#         cluster="mapresults/Genecatalog_CE/canopy_cluster.tsv",
#         profile="mapresults/Genecatalog_CE/cluster_profiles.tsv",
#     params:
#         canopy_params=config.get("canopy_params","")
#     log:
#         "mapresults/Genecatalog_CE/canopy.log"
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
