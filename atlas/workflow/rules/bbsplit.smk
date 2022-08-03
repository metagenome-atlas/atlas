


def run_bbsplit(*args,log=None,**kwargs):

    command='bbsplit.sh'

    for a in args:
        command+=' '+a
    for k in kwargs:
        if type(kwargs[k])==bool:
            kwargs[k]= 't' if kwargs[k] else 'f'
        command+=f" {k}={kwargs[k]} "
    if log is not None:
        command += f' 2>{log}'

    logger.info(f"run: {command}")
    shell(command)



rule bbsplit_index:
    input:
        genome_dir,
    output:
        index= directory("genomes/bbsplit_ref/ref"),
    threads:
        config.get("threads", 6)
    resources:
        mem_mb = config["mem"] * 1000,
        java_mem = int(config["mem"] * JAVA_MEM_FRACTION),
        time_min = 60 * config['runtime']["long"]
    log:
        "logs/genomes/mapping/build_bbmap_index.log",
    benchmark:
        "logs/benchmark/genomes/bbsplit_index.log"
    run:
        run_bbsplit(f"-Xmx{resources.java_mem}G" ,
                    ref= input,
                    path= os.path.dirname(output[0]),
                    threads=threads,
                    log=log[0]
                    )




# generalized rule so that reads from any "sample" can be aligned to contigs from "sample_contigs"
rule bbsplit:
    input:
        reads= lambda wildcards: get_quality_controlled_reads_(wildcards,['R1','R2']),
        refdir = rules.bbsplit_index.output.index,
    output:
        sam=temp("genomes/alignments/sam/{sample}.sam"),
        scafstats = temp("genomes/alignments/scafstats/{sample}.tsv.gz"),
        refstats = temp("genomes/alignments/refstats/{sample}.tsv.gz"),
        covstats = "genomes/alignments/coverage/{sample}.tsv.gz",
        bincov = temp("genomes/alignments/coverage_binned/{sample}.tsv.gz")
    params:
        ambiguous="all" ,
        ambiguous2 = "best",
        minid = config.get('contig_min_id', CONTIG_MIN_ID),
        #unmapped=lambda wc, output: io_params_for_tadpole(output.unmapped, "outu"),
    shadow:
        "shallow"        
    log:
        "logs/genomes/bbsplit/{sample}.log"
    benchmark:
        "logs/benchmark/genomes/bbsplit/{sample}.txt"
    conda:
        "../envs/required_packages.yaml"
    threads:
        config["threads"]
    resources:
        mem_mb = 1000 * config.get("java_mem", JAVA_MEM),
        java_mem = int(config.get("java_mem", JAVA_MEM) * JAVA_MEM_FRACTION),
        time_min= config["runtime"]["default"] * 60
        time= config["runtime"]["default"]
    run:

        run_bbsplit(f"-Xmx{resources.java_mem}G" ,
                    io_params_for_tadpole(input.reads,'in'),
                    path= os.path.dirname(input.refdir),
                    threads=threads,
                    log=log[0],
                    **output,
                    **params
                    )


localrules: count_mapped_reads



rule count_mapped_reads:
    input:
        logfile=expand( "logs/genomes/bbsplit/{sample}.log", sample=SAMPLES),
        out_check=expand("genomes/alignments/scafstats/{sample}.tsv.gz", sample=SAMPLES)
    output:
        "genomes/counts/mapping_rate.tsv"
    run:
        import pandas as pd
        from utils.parsers_bbmap import parse_bbmap_log_file

        D= pd.DataFrame(index=SAMPLES,columns=['reads_used','reads_mapped'])
        for i,sample in enumerate(SAMPLES):
            D.loc[sample]= parse_bbmap_log_file(input[i])

        D['mapping_rate']= D.iloc[:,1]/D.iloc[:,0]
        D.to_csv(output[0],sep='\t')


rule combine_bined_coverages_MAGs:
    input:
        binned_coverage_files=expand(
           "genomes/alignments/coverage_binned/{sample}.tsv.gz", sample=SAMPLES
        ),
        
    params:
        samples=SAMPLES,
    output:
        binned_cov="genomes/counts/binned_coverage.parquet",
        median_abund="genomes/counts/median_coverage_genomes.parquet",
    log:
        "logs/genomes/counts/combine_binned_coverages_MAGs.log",
    threads: 1
    script:
        "../scripts/combine_coverage_MAGs.py"


rule merge_counts:
    input:
        expand("genomes/alignments/{{scope}}stats/{sample}.tsv.gz",sample=SAMPLES)
    output:
        counts="genomes/counts/combined_coverage_{scope}.parquet",
    run:
        import pandas as pd
        import gc
        samples= SAMPLES

        Reads={}
        for i in range(len(input)):
            df=pd.read_csv(input[i],index_col=0,sep='\t')
            Reads[samples[i]]= pd.to_numeric(df.unambiguousReads,downcast='usigned')
            del df
            gc.collect()


        Reads= pd.concat(Reads,axis=1,sort=False).fillna(0).T
        gc.collect()

        Reads.reset_index().to_parquet(output.counts)
        
