import os



if config['genecatalog']['source']=='contigs':

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

else:

    localrules: concat_genes
    rule concat_genes:
        input:
            "genomes/annotations/genes"
        output:
            faa=  temp("Genecatalog/all_genes_unfiltered.faa"),
            fna = temp("Genecatalog/all_genes_unfiltered.fna"),
        shell:
            " cat {input}/*.faa >  {output.faa} ;"
            " cat {input}/*.fna > {output.fna}"


localrules: filter_genes
rule filter_genes:
    input:
        fna="Genecatalog/all_genes_unfiltered.fna",
        faa="Genecatalog/all_genes_unfiltered.faa"
    output:
        fna= "Genecatalog/all_genes/predicted_genes.fna",
        faa= "Genecatalog/all_genes/predicted_genes.faa",
    threads:
        1
    params:
        min_length=config['genecatalog']['minlength']
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


if (config['genecatalog']['clustermethod']=='linclust') or (config['genecatalog']['clustermethod']=='mmseqs'):

    rule cluster_genes:
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
            tmpdir= os.path.join(config['tmpdir'],"mmseqs"),
            clustermethod = 'linclust' if config['genecatalog']['clustermethod']=='linclust' else 'cluster',
            coverage=config['genecatalog']['coverage'], #0.8,
            minid=config['genecatalog']['minid'], # 0.00
            extra=config['genecatalog']['extra']
        shell:
            """
                mmseqs createdb {input.faa} {output.db[0]} > >(tee  {log})

                mkdir -p {params.tmpdir}

                mmseqs {params.clustermethod} -c {params.coverage} \
                --min-seq-id {params.minid} {params.extra} \
                --threads {threads} {output.db[0]} {output.clusterdb[0]} {params.tmpdir}  > >(tee -a  {log})

                rm -fr  {params.tmpdir} > >(tee -a  {log})
            """


    rule get_rep_proteins:
        input:
            db= rules.cluster_genes.output.db,
            clusterdb = rules.cluster_genes.output.clusterdb,
        output:
            cluster_attribution = temp("Genecatalog/orf2gene_oldnames.tsv"),
            rep_seqs_db = temp(expand("Genecatalog/protein_catalog.{exp}",exp=['db','db.index'])),
            rep_seqs = temp("Genecatalog/representatives_of_clusters.fasta")
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


    localrules: rename_protein_catalog
    rule rename_protein_catalog:
        input:
            cluster_attribution = "Genecatalog/orf2gene_oldnames.tsv",
        output:
            cluster_attribution = "Genecatalog/clustering/orf2gene.tsv",
        run:
            import pandas as pd
            # CLuterID    GeneID    empty third column
            gene2proteins= pd.read_table(input.cluster_attribution,index_col=1, header=None)

            protein_clusters_old_names= gene2proteins[0].unique()

            map_names = dict(zip(protein_clusters_old_names,
                                 gen_names_for_range(len(protein_clusters_old_names),'Gene')))

            gene2proteins['Gene'] = gene2proteins[0].map(map_names)
            gene2proteins.index.name='ORF'
            gene2proteins['Gene'].to_csv(output.cluster_attribution,sep='\t',header=True)



elif config['genecatalog']['clustermethod']=='cd-hit-est':

# cluster genes


    rule cluster_genes:
        input:
            fna_dir="Genecatalog/all_genes/predicted_genes.fna",
        output:
            temp("Genecatalog/representatives_of_clusters.fasta"),
            temp("Genecatalog/gene_catalog_oldnames.clstr")
        conda:
            "%s/cd-hit.yaml" % CONDAENV
        log:
            "logs/Genecatalog/cluster_genes.log"
        threads:
            config.get("threads", 1)
        resources:
            mem= config.get("java_mem", JAVA_MEM)
        params:
            coverage=config['genecatalog']['coverage'],
            identity=config['genecatalog']['minid'],
            extra= config['genecatalog']['extra'],
            prefix= lambda wc,output: os.path.splitext(output[1])[0],
        shell:
            """
                cd-hit-est -i {input} -T {threads} \
                -M {resources.mem}000 -o {params.prefix} \
                -c {params.identity} -n 9  -d 0 {params.extra} \
                -aS {params.coverage} -aL {params.coverage} > >(tee {log})

                mv {params.prefix} {output[0]} 2>> {log}
            """





    localrules:  parse_clstr_files, rename_gene_clusters

    def parse_cd_hit_file(clstr_file):
        """

        >Cluster 0
        0	342nt, >S1_83_1... *
        1	342nt, >S2_82_1... at +/100.00%
        >Cluster 1
        0	339nt, >S1_61_1... *
        1	339nt, >S2_59_1... at +/100.00%


        """
        import numpy as np
        def parse_line(line):
            _, length, name, identity = line.strip().replace('...','\t').replace(', ','\t').split('\t')

            length= int(length.replace('nt',''))
            name=name[1:]
            if '*' in identity:
                identity= np.nan
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
                    if np.isnan(identity):
                        cluster['representative']= name
        return Clusters


    def write_cd_hit_clusters(Clusters,file_handle):
            for cluster in Clusters:
                for element in cluster['elements']:
                    file_handle.write(f"{element[0]}\t{element[1]}\t{element[2]}\t{cluster['representative']}\n")



    rule parse_clstr_files:
        input:
            clustered_dir= "Genecatalog/gene_catalog_oldnames.clstr"
        output:
            temp("Genecatalog/clustering/orf2gene_oldnames.tsv")
        run:
            with open(output[0],'w') as fout:
                fout.write(f"ORF\tLength\tIdentity\tGene\n")
                Clusters=parse_cd_hit_file(input[0])
                write_cd_hit_clusters(Clusters,fout)



    rule rename_gene_clusters:
        input:
            orf2gene = "Genecatalog/clustering/orf2gene_oldnames.tsv",
        output:
            orf2gene = "Genecatalog/clustering/orf2gene.tsv",
        run:
            import pandas as pd
            from Bio import SeqIO

            orf2gene= pd.read_table(input.orf2gene,index_col=0)

            # rename gene repr to Gene0000XX

            gene_clusters_old_names= orf2gene['Gene'].unique()

            map_names = dict(zip(gene_clusters_old_names,
                                 gen_names_for_range(len(gene_clusters_old_names),'Gene')))

            orf2gene['Gene'] = orf2gene['Gene'].map(map_names)
            orf2gene.to_csv(output.orf2gene,sep='\t',header=True)

else:
    raise Exception("Didn't understood the genecatalog clustermethod: {}".format(config['genecatalog']['clustermethod']))

localrules: rename_gene_catalog
rule rename_gene_catalog:
    input:
        fna = "Genecatalog/all_genes/predicted_genes.fna",
        faa= "Genecatalog/all_genes/predicted_genes.faa",
        orf2gene = "Genecatalog/clustering/orf2gene.tsv",
        representatives= "Genecatalog/representatives_of_clusters.fasta"
    output:
        fna= "Genecatalog/gene_catalog.fna",
        faa= "Genecatalog/gene_catalog.faa",
    run:
        import pandas as pd
        from Bio import SeqIO

        representatives= []
        with open(input.representatives) as fasta:
            for line in fasta:
                if line[0]=='>': representatives.append(line[1:].split()[0])

        map_names= pd.read_table(input.orf2gene,index_col=0).loc[representatives,'Gene']

        # rename fna
        faa_parser = SeqIO.parse(input.faa,'fasta')
        fna_parser = SeqIO.parse(input.fna,'fasta')

        with open(output.fna,'w') as fna, open(output.faa,'w') as faa :
            for gene in fna_parser:
                protein = next(faa_parser)
                if gene.name in map_names.index:
                    gene.id = map_names[gene.name]
                    protein.id = map_names[protein.name]

                    SeqIO.write(gene,fna,'fasta')
                    SeqIO.write(protein,faa,'fasta')


rule align_reads_to_Genecatalog:
    input:
        unpack(get_quality_controlled_reads),
        fasta = "Genecatalog/gene_catalog.fna",
    output:
        sam = temp("Genecatalog/alignments/{sample}.sam")
    params:
        input = lambda wc, input : input_params_for_bbwrap(wc, input),
        maxsites = 4,
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
            out={output.sam} \
            threads={threads} \
            minid={params.minid} \
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
        "Genecatalog/counts/median_coverage.tsv.gz",
        "Genecatalog/counts/Nmapped_reads.tsv.gz",
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

        pd.DataFrame(combined_cov).to_csv(output[0],sep='\t',compression='gzip')
        pd.DataFrame(combined_N_reads).to_csv(output[1],sep='\t',compression='gzip')


###########
## EGG NOG
##########

# # this rule specifies the more general eggNOG rules

# output with wildcards "{folder}/{prefix}.emapper.tsv"


def get_eggnog_db_file():
    return expand("{path}/{files}",
                  path=EGGNOG_DIR,
                  files=["OG_fasta","eggnog.db","og2level.tsv","eggnog_proteins.dmnd"]
                  )

# TODO: make benchmark
rule eggNOG_homology_search:
    input:
        eggnog_db_files=get_eggnog_db_file(),
        faa = "{folder}/{prefix}.faa",
    output:
        temp("{folder}/{prefix}.emapper.seed_orthologs"),
    params:
        data_dir = EGGNOG_DIR,
        prefix = "{folder}/{prefix}"
    resources:
        mem = config.get("java_mem", JAVA_MEM)
    threads:
        config["threads"]
    conda:
        "%s/eggNOG.yaml" % CONDAENV
    log:
        "{folder}/logs/{prefix}/eggNOG_homology_search_diamond.log"
    shell:
        """
        emapper.py -m diamond --no_annot --no_file_comments \
            --data_dir {params.data_dir} --cpu {threads} -i {input.faa} \
            -o {params.prefix} --override 2> >(tee {log})
        """



rule eggNOG_annotation:
    input:
        eggnog_db_files=get_eggnog_db_file(),
        seed = rules.eggNOG_homology_search.output
    output:
        temp("{folder}/{prefix}.emapper.annotations")
    params:
        data_dir = EGGNOG_DIR,
        prefix = "{folder}/{prefix}"
    threads:
        config.get("threads", 1)
    resources:
        mem=20
    conda:
        "%s/eggNOG.yaml" % CONDAENV
    log:
        "{folder}/logs/{prefix}/eggNOG_annotate_hits_table.log"
    shell:
        """
        emapper.py --annotate_hits_table {input.seed} --no_file_comments --usemem \
            --override -o {params.prefix} --cpu {threads} --data_dir {params.data_dir} 2> >(tee {log})
        """

EGGNOG_HEADERS= [
"query_name",
"seed_eggNOG_ortholog",
"seed_ortholog_evalue",
"seed_ortholog_score",
"predicted_gene_name",
"GO_terms",
"KEGG_KO",
"BiGG_Reactions",
"Annotation_tax_scope",
"Matching_OGs",
"best_OG|evalue|score",
"categories",
"eggNOG_HMM_model_annotation"]

# rule add_eggNOG_header:
#     input:
#         "{folder}/{prefix}.emapper.annotations"
#     output:
#         "{folder}/{prefix}.emapper.tsv"
#     run:
#         import pandas as pd#

#            where do you take the Headers

#         D = pd.read_table(input[0], header=None)
#         D.columns = EGGNOG_HEADERS
#         D.to_csv(output[0],sep="\t",index=False)





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
        "Genecatalog/gene_catalog.faa"
    output:
        "Genecatalog/annotation/single_copy_genes_{domain}.tsv",
    params:
        script_dir = os.path.dirname(os.path.abspath(workflow.snakefile)),
        key = lambda wc: wc.domain[:3] #bac for bacteria, #arc for archaea
    conda:
        "%s/DASTool.yaml" % CONDAENV # needs pearl
    log:
        "logs/Genecatalog/annotation/predict_single_copy_genes_{domain}.log"
    shadow:
        "shallow"
    threads:
        config['threads']
    shell:
        " DIR=$(dirname $(readlink -f $(which DAS_Tool))) "
        ";"
        " ruby {params.script_dir}/rules/scg_blank_diamond.rb diamond"
        " {input} "
        " $DIR\/db/{params.key}.all.faa "
        " $DIR\/db/{params.key}.scg.faa "
        " $DIR\/db/{params.key}.scg.lookup "
        " {threads} "
        " 2> >(tee {log}) "
        " ; "
        " mv {input[0]}.scg {output}"




localrules: gene_subsets,combine_egg_nogg_annotations
checkpoint gene_subsets:
    input:
        "Genecatalog/gene_catalog.faa"
    output:
        directory("Genecatalog/subsets/genes")
    params:
        subset_size=config['genecatalog']['SubsetSize'],
    run:
        from utils import fasta
        fasta.split(input[0],params.subset_size,output[0],simplify_headers=True)


def combine_genecatalog_annotations_input(wildcards):
    dir_for_subsets = checkpoints.gene_subsets.get(**wildcards).output[0]
    Subset_names,= glob_wildcards(os.path.join(dir_for_subsets, "{subset}.faa"))
    return expand("Genecatalog/subsets/genes/{subset}.emapper.annotations",
                  subset=Subset_names)

rule combine_egg_nogg_annotations:
    input:
        combine_genecatalog_annotations_input
    output:
        temp("Genecatalog/annotations/eggNog.emapper.annotations")
    shell:
        "cat {input} > {output}"

localrules: add_eggNOG_header
rule add_eggNOG_header:
    input:
        "Genecatalog/annotations/eggNog.emapper.annotations"
    output:
        "Genecatalog/annotations/eggNog.tsv"
    run:
        import pandas as pd

        D = pd.read_table(input[0], header=None)
        D.columns = EGGNOG_HEADERS
        D.to_csv(output[0],sep="\t",index=False)

# localrules: gene_subsets
# rule gene_subsets:
#     input:
#         "Genecatalog/gene_catalog.faa"
#     output:
#         dynamic("Genecatalog/subsets/genes/{subset}.faa")
#     params:
#         subset_size=config['genecatalog']['SubsetSize'],
#     run:
#         from utils import fasta
#
#         output_dir=os.path.dirname(output[0])
#         os.removedirs(output_dir)
#
#         fasta.split(input[0],params.subset_size,output_dir,simplify_headers=True)
#
# rule combine_annotations:
#     input:
#         eggNOG=dynamic("Genecatalog/subsets/genes/{subset}.emapper.annotations")
#     output:
#         eggNOG= "Genecatalog/annotations/eggNog.tsv"
#     run:
#
#         with open(input[0]) as f:
#             first_line= f.readline()
#             assert len(first_line.split('\t')) == len(EGGNOG_HEADERS), "number of eggnog headers doesn't correspond to number of fields."
#
#         with open(output.eggNOG,'w') as f:
#             f.write("\t".join(EGGNOG_HEADERS) + '\n')
#         shell("cat {input.eggNOG} >> {output.eggNOG}")
#




# after combination need to add eggNOG headerself.
#"{folder}/{prefix}_eggNOG.tsv"



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
