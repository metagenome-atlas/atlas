


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
