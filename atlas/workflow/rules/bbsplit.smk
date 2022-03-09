
import os
import re
import sys
from glob import glob

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)),"scripts"))


import warnings

include: "sample_table.smk"

MAXIMUM_COUNTED_MAP_SITES=10
CONTIG_MAP_PAIRED_ONLY = True
CONTIG_COUNT_MULTI_MAPPED_READS=False
CONTIG_MAX_DISTANCE_BETWEEN_PAIRS = 300
CONTIG_MIN_ID= 0.9
JAVA_MEM=250
config['java_mem']=JAVA_MEM
config['index_mem']=250
config['index_time']=12
config['threads']= 8
config['ambiguous2']='best'

if not ('ref_name' in config):
    config['ref_name']='uhgg2'

config['ref_paths']={'strains': '/home/kiesers/scratch/MAGcatalog/representatives/strains' ,
                     'species': '/home/kiesers/scratch/MAGcatalog/representatives/species',
                     #'uhgg2': '/home/kiesers/scratch/HumanUnifiedCatalog/representatives/unzip/',
                     'mMAG': "/home/kiesers/scratch/iMGMC/dereplicated_genomes"
                    }

JAVA_MEM_FRACTION=0.85
CONDAENV='envs'

Paired_Single=[]
if PAIRED_END:
    Paired_Single.append('pe')
else:
#if 'se' in MULTIFILE_FRACTIONS:
    Paired_Single.append('se')

print(f"Use reference {config['ref_name']}")

def get_reference_input(wildcards):

    ref_path= os.path.abspath(config['ref_paths'][wildcards.ref_name])
    assert os.path.exists(ref_path), f"Reference doesn't exist: {ref_path}"
    assert os.path.isdir(ref_path), f"Reference should be a dir with fasta files {ref_path}"

    return ref_path



# def get_index_path(wildcards):
#     return os.path.join( get_reference_input(wildcards),"../bbsplit_ref",wildcards.ref_name,"ref")


DBDIR = config.get("database_dir",'/home/kiesers/scratch/Atlas/databases/bbsplit_references')

logger.info(f"{len(SAMPLES)} samples")

wildcard_constraints:
    sample="\w+",
    ref_name="\w+"


rule all:
    input:
        f"mapping_results/{config['ref_name']}/mapping_rate.tsv",
        expand("mapping_results/{ref_name}/counts/median_coverage_genomes_{pese}.tsv.gz",pese=Paired_Single,ref_name=config['ref_name']) ,
        expand("mapping_results/{ref_name}/counts/counts_{fraction}_{unit}.tsv.gz",
               unit=['ref','scaf'], fraction=Paired_Single,ref_name=config['ref_name']),
        #f"{os.path.abspath( config['ref_paths'][ config['ref_name'] ] )}/../bbsplit_ref/contigs2{config['ref_name']}.tsv"




def run_bbsplit(*args,log=None,**kwargs):

    command='bbsplit.sh'

    for a in args:
        command+=' '+a
    for k in kwargs:
        if type(kwargs[k])==bool:
            kwargs[k]= 't' if kwargs[k] else 'f'
        command+=f" {k}={kwargs[k]} "
    if log is not None:
        if kwargs.get('append',False):
            command += f' 2>>{log}'
        else:
            command += f' 2>{log}'

    logger.info(f"run: {command}")
    shell(command)



rule bbsplit_index:
    input:
        get_reference_input,
    output:
        index= directory(f"{DBDIR}/{{ref_name}}/ref"),
    threads:
        config.get("threads", 6)
    resources:
        mem_mb = 1000 * config["index_mem"],
        java_mem = int(config["index_mem"] * JAVA_MEM_FRACTION),
        time_min = 60 * config['index_time']
    log:
        f"{DBDIR}/{{ref_name}}/bbsplit_index.log"
    benchmark:
        f"{DBDIR}/{{ref_name}}/bbsplit_index_benchmark.log"
    run:
        run_bbsplit(f"-Xmx{resources.java_mem}G" ,
                    ref= input,
                    path= os.path.dirname(output[0]),
                    threads=threads,
                    log=log[0]
                    )




# generalized rule so that reads from any "sample" can be aligned to contigs from "sample_contigs"
rule bbsplit_pe:
    input:
        reads= lambda wildcards: get_quality_controlled_reads_(wildcards,['R1','R2']),
        refdir = f"{DBDIR}/{{ref_name}}/ref",
    output:
        #sam = temp("{sample}.sam"),
        #unmapped = expand("mapping_results/{{ref_name}}/unmapped/{{sample}}_{fraction}.fastq.gz",
        #                  fraction=['R1','R2']),
        scafstats="mapping_results/{ref_name}/counts/pe/{sample}_scafstats.tsv.gz",
        refstats= "mapping_results/{ref_name}/counts/pe/{sample}_refstats.tsv.gz",
        covstats = "mapping_results/{ref_name}/counts/pe/{sample}_coverage.tsv.gz",
        bincov = "mapping_results/{ref_name}/counts/pe/{sample}_coverage_binned.tsv.gz"
    params:
        ambiguous2 = config['ambiguous2'],
        minid = config.get('contig_min_id', CONTIG_MIN_ID),
    log:
        "logs/bbsplit/{ref_name}/pe/{sample}.log"
    benchmark:
        "logs/benchmark/bbsplit/{ref_name}/{sample}_pe.txt"
#    conda:
#        "%s/required_packages.yaml" % CONDAENV
    threads:
        config["threads"]
    resources:
        mem_mb = 1000 * config.get("java_mem", JAVA_MEM),
        java_mem = int(config.get("java_mem", JAVA_MEM) * JAVA_MEM_FRACTION),
        time_min=12*60
    run:

        run_bbsplit(f"-Xmx{resources.java_mem}G" ,
                    io_params_for_tadpole(input.reads,'in'),
                    path= os.path.dirname(input.refdir),
                    threads=threads,
                    log=log[0],
                    **output,
                    **params
                    )

rule bbsplit_se:
    input:
        reads= lambda wildcards: get_quality_controlled_reads_(wildcards,['se']),
        refdir = f"{DBDIR}/{{ref_name}}/ref",
    output:
        #sam = temp("{sample}.sam"),
        #unmapped = "mapping_results/{ref_name}/unmapped/{sample}_se.fastq.gz",
        scafstats=temp("mapping_results/{ref_name}/counts/se/{sample}_scafstats.tsv.gz"),
        refstats= temp("mapping_results/{ref_name}/counts/se/{sample}_refstats.tsv.gz"),
        covstats = temp("mapping_results/{ref_name}/counts/se/{sample}_coverage.tsv.gz"),
        bincov = temp("mapping_results/{ref_name}/counts/se/{sample}_coverage_binned.tsv.gz")
    params:
        ambiguous2 = config['ambiguous2'],
        minid = config.get('contig_min_id', CONTIG_MIN_ID),
        physcov=True
    log:
        "logs/bbsplit/{ref_name}/se/{sample}.log"
    benchmark:
        "logs/benchmark/bbsplit/{ref_name}/{sample}_se.txt"
#    conda:
#        "%s/required_packages.yaml" % CONDAENV
    threads:
        config["threads"]
    resources:
        mem_mb = 1000* config.get("java_mem", JAVA_MEM),
        java_mem = int(config.get("java_mem", JAVA_MEM) * JAVA_MEM_FRACTION),
        time_min=8*60
    run:

        run_bbsplit(f"-Xmx{resources.java_mem}G" ,
                    io_params_for_tadpole(input.reads,'in'),
                    path= os.path.dirname(input.refdir),
                    threads=threads,
                    log=log[0],
                    **output,
                    **params
                    )




localrules: count_mapped_reads, merge_counts,merge_log_files

rule merge_log_files:
    input:
        log=expand("logs/bbsplit/{{ref_name}}/{pese}/{{sample}}.log", pese=Paired_Single),
        out_check=expand("mapping_results/{{ref_name}}/counts/{pese}/{{sample}}_refstats.tsv.gz", pese=Paired_Single)
    output:
        "logs/bbsplit/{ref_name}/{sample}.log"
    shell:
        "cat {input.log} > {output}"

#only use paired end reads

rule count_mapped_reads:
    input:
        #logfile=expand("logs/bbsplit/{{ref_name}}/{sample}.log",sample=SAMPLES),
        logfile=expand("logs/bbsplit/{{ref_name}}/{pese}/{sample}.log", pese=Paired_Single,sample=SAMPLES),
        out_check=expand("mapping_results/{{ref_name}}/counts/{pese}/{sample}_refstats.tsv.gz", pese=Paired_Single,sample=SAMPLES)
    output:
        "mapping_results/{ref_name}/mapping_rate.tsv"
    run:
        import pandas as pd
        from parsers_bbmap import parse_bbmap_log_file

        D= pd.DataFrame(index=SAMPLES,columns=['reads_used','reads_mapped'])
        for i,sample in enumerate(SAMPLES):
            D.loc[sample]= parse_bbmap_log_file(input[i])

        D['mapping_rate']= D.iloc[:,1]/D.iloc[:,0]
        D.to_csv(output[0],sep='\t')



rule merge_counts:
    input:
        expand("mapping_results/{{ref_name}}/counts/{{pese}}/{sample}_{{aggregation}}stats.tsv.gz",sample=SAMPLES)
    output:
        counts="mapping_results/{ref_name}/counts/counts_{pese}_{aggregation}.tsv.gz",
        #relab= "mapping_results/{ref_name}/counts/relab_{pese}_{aggregation}.tsv.gz"
    run:
        import pandas as pd

        samples= SAMPLES

        Reads={}
        #Relab= {}
        for i in range(len(input)):
            df=pd.read_csv(input[i],index_col=0,sep='\t')
            Reads[samples[i]]= df.unambiguousReads
            #Relab[samples[i]]= df['%unambiguousReads']/100

        Reads= pd.concat(Reads,axis=1,sort=False).fillna(0).T
        #Relab= pd.concat(Relab,axis=1,sort=False).fillna(0).T

        Reads.to_csv(output.counts,sep='\t')
        #Relab.to_csv(output.relab)



# rule combine_coverages:
#     input:
#         covstats = expand("mapping_results/{{ref_name}}/counts/{{pese}}/{sample}_coverage.txt.gz",
#             sample=SAMPLES),
#     output:
#         "mapping_results/{ref_name}/counts/median_contig_coverage_{pese}.tsv.gz",
#         "mapping_results/{ref_name}/counts/raw_counts_contigs_{pese}.tsv.gz"
#     run:
#         import pandas as pd
#         from parsers_bbmap import combine_coverages
#
#         combined_cov,Counts_contigs = combine_coverages(input.covstats,SAMPLES)
#
#         combined_cov.to_csv(output[0],sep='\t')
#         Counts_contigs.to_csv(output[1],sep='\t')





# [identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment]
localrules: get_contigs2bin
rule get_contigs2bin:
    input:
        "{referece_path}"
    output:
        "{referece_path}/../bbsplit_ref/contigs2{ref_name}.tsv",
    run:

        import os
        from glob import glob
        mapfile_contigs= output[0]
        fasta_files= glob(input[0]+'/*.f*')

        with open(mapfile_contigs,'w') as out_contigs :
            for fasta in fasta_files:

                bin_name= os.path.splitext(os.path.split(fasta)[-1])[0]

                # write names of contigs in mapping file
                with open(fasta) as f :
                    for line in f:
                        if line[0]==">":
                            header=line[1:].strip()
                            out_contigs.write(f"{header}\t{bin_name}\n")

rule combine_bined_coverages:
    input:
        binned_coverage_files = expand("mapping_results/{{ref_name}}/counts/{{pese}}/{sample}_coverage_binned.tsv.gz",
            sample=SAMPLES),
    params:
        samples= SAMPLES
    output:
        median_abund = "mapping_results/{ref_name}/counts/median_coverage_genomes_{pese}.tsv.gz"
    resources:
        mem_mb = 250 * 1000,
        time_min=8 *60
    benchmark:
        "logs/benchmark/bbsplit/{ref_name}/combined_binned_coverages_{pese}.txt"
    run:

        import pandas as pd
        import os


        def read_median_abundance(covarage_binned_file):
            """
            #Mean   1.545
            #STDev  29.298
            #RefName        Cov     Pos     RunningPos
            MAG20067$ERR675668_24   0.00    1000    0
            MAG20067$ERR675668_24   0.00    2000    1000

            """
            binCov= pd.read_csv(covarage_binned_file,sep='\t',
                             skiprows=2,
                             index_col=[0,2],
                             usecols=[0,1,2],
                             squeeze=True
                             )
            # split first index in two genome$contig

            index= pd.Series(binCov.index.levels[0],index=binCov.index.levels[0] )
            splitted= index.str.split('$',expand=True)
            splitted.columns=['Genome','Contig']
            new_index= splitted.loc[binCov.index.get_level_values(0)]
            new_index['Position']= binCov.index.get_level_values(1).values
            binCov.index=pd.MultiIndex.from_frame(new_index)

            median_abund= binCov.groupby(level=0).median().T

            return median_abund



        Median_abund={}
        for i, cov_file in enumerate(input.binned_coverage_files):

            sample= params.samples[i]

            print(f'read sample: {sample}')
            Median_abund[sample] = read_median_abundance(cov_file)

        print('combine samples')
        Median_abund = pd.concat(Median_abund)

        Median_abund.to_csv(output.median_abund,sep='\t')


# rule combine_bined_coverages:
#     input:
#         binned_coverage_files = expand("mapping_results/{{ref_name}}/counts/{{pese}}/{sample}_coverage_binned.tsv.gz",
#             sample=SAMPLES),
#     params:
#         samples= SAMPLES
#     output:
#         binned_cov= "mapping_results/{ref_name}/counts/combined_binned_coverage_{pese}.tsv.gz",
#         median_abund = "mapping_results/{ref_name}/counts/median_coverage_genomes_{pese}.tsv.gz"
#     resources:
#         mem_mb = 250*1000,
#         time_min=8*60
#     benchmark:
#         "logs/benchmark/bbsplit/{ref_name}/combined_binned_coverages_{pese}.txt"
#     run:
#
#         import pandas as pd
#         import os
#
#
#         def read_coverage_binned(covarage_binned_file):
#             """
#             #Mean   1.545
#             #STDev  29.298
#             #RefName        Cov     Pos     RunningPos
#             MAG20067$ERR675668_24   0.00    1000    0
#             MAG20067$ERR675668_24   0.00    2000    1000
#
#             """
#             return pd.read_csv(covarage_binned_file,sep='\t',
#                              skiprows=2,
#                              index_col=[0,2],
#                              usecols=[0,1,2],
#                              squeeze=True)
#
#
#
#         binCov={}
#         for i, cov_file in enumerate(input.binned_coverage_files):
#
#             sample= params.samples[i]
#
#             print(f'read sample: {sample}')
#             binCov[sample] = read_coverage_binned(cov_file)
#
#         print('combine samples')
#         binCov = pd.concat(binCov,copy=False)
#
#         # split first index in two genome$contig
#
#         index= pd.Series(binCov.index.levels[0],index=binCov.index.levels[0] )
#         splitted= index.str.split('$',expand=True)
#         splitted.columns=['Genome','Contig']
#         new_index= splitted.loc[binCov.index.get_level_values(0)]
#         new_index['Position']= binCov.index.get_level_values(1).values
#         binCov.index=pd.MultiIndex.from_frame(new_index)
#
#         print('save combined table')
#         binCov.to_csv(output.binned_cov,sep='\t',compression='gzip')
#
#         print('Calculate median')
#         Median_abund= binCov.groupby(level=0).median().T
#
#         Median_abund.to_csv(output.median_abund,sep='\t')
