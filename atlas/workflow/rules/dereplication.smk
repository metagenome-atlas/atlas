


rule bindash_sketch_genome:
    input:
        "genomes/all_bins/filtered_bins.txt"
    output:
        "genomes/sketches/genomes.bdsh",
        "genomes/sketches/genomes.bdsh.dat",
        "genomes/sketches/genomes.bdsh.txt"
    params:
        k= config['sketch_k'],
        sketchsize64= int(config['sketch_size'])//64,
        extra=config.get('bindash_extra',"")
    threads:
        config['threads']
    conda:
        "../envs/bindash.yaml"
    log:
        "logs/bindash/sketch.log"
    benchmark:
        "logs/benchmark/bindash_sketch.txt"
    shell:
        "bindash sketch "
        "--outfname={output[0]} "
        "--nthreads={threads} "
        "--sketchsize64={params.sketchsize64} "
        "--kmerlen={params.k} "
        "{params.extra} "
        "--listfname={input[0]} 2> {log}"

rule bindash_dist:
    input:
        rules.bindash_sketch_genome.output[0]
    output:
        "genomes/sketches/bindash_dists.tsv"
    params:
        d= config['sketch_max_dist']
    threads:
        config['threads']
    conda:
        "../envs/bindash.yaml"
    log:
        "logs/bindash/dist.log"
    benchmark:
        "logs/benchmark/bindash_dist.txt"
    shell:
        "bindash dist "
        "--nthreads={threads} "
        "--mthres={params.d} "
        "--outfname={output} {input[0]} 2> {log}"

rule tsv2parquet:
    input:
        rules.bindash_dist.output
    output:
        "genomes/clustering/sketch_dists.parquet"
    resources:
        mem_mb=config['mem']['large'] *1000
    threads:
        1
    run:

        from utils.sketches import load_bindash

        M= open_function(input[0]).drop(['Identity'],axis=1)
        M.to_parquet(output[0],engine="pyarrow")

# TODO: adapt this rule
checkpoint cluster_species:
    input:
        dists= rules.tsv2parquet.output,
        genome_info = "genomes/all_bins/filtered_quality.tsv",
    output:
        cluster_file="genomes/clustering/bins_clustering.tsv",
    resources:
        mem_mb=config['mem']['large']*1000
    params:
        threshold=config["genome_dereplication"]["ANI"],
        linkage_method= 'average',
    script:
        "../scripts/cluster_genomes.py"


def get_species(wildcards):
    import pandas as pd
    cluster_file=checkpoints.cluster_species.get().output.cluster_file

    df= pd.read_csv(cluster_file,sep='\t',index_col=0)
    return list(df.Species.unique())