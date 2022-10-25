


rule bindash_sketch_genome:
    input:
        "sketches/genome_list_bindash.txt"
    output:
        "sketches/genomes.bdsh",
        "sketches/genomes.bdsh.dat",
        "sketches/genomes.bdsh.txt"
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
        "tables/bindash_dists.tsv"
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
        "tables/{tool}_dists.tsv"
    output:
        "tables/{tool}_dists.parquet"
    resources:
        mem_mb=config['mem']['large'] *1000
    threads:
        1
    run:

        if wildcards.tool == "mummer":
            open_function= gd.load_mummer
        elif wildcards.tool == "minimap":
            open_function= gd.load_minimap
        elif wildcards.tool == "bindash":
            open_function= gd.load_bindash
        elif wildcards.tool == "fastani":
            open_function= gd.load_fastani
        elif wildcards.tool == "mash":
            open_function= gd.load_mash
        else:
            raise Exception(
                f"Don't know how to load table from tool : {wildcards.tool}"
            )


        M= open_function(input[0]).drop(['Identity'],axis=1)
        M.to_parquet(output[0],engine="pyarrow")
