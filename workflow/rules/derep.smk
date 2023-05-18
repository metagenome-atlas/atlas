


rule run_skani:
    input:
        paths="Binning/{binner}/filtered_bins_paths.txt"
    output:
        temp("Intermediate/dereplication/{binner}_distance_matrix.txt")
    log:
        "logs/binning/dereplication/skani_calculation_{binner}.log"
    resources:
        mem_mb= config["mem"]*1000 ,
        time_min=60 * config["runtime"]["default"],
    params:
        preset= "medium" # fast, medium or slow
    threads:
        config["threads"],
    conda:
        "../envs/skani.yaml"
    shell:
        "skani triangle "
        " -l {input.paths} "
        " -o {output} "
        " -t {threads} "
        " --trace "
        " --robust "
        " --sparse --ci "
        " --min-af 15"
        " &> {log} "



rule skani_2_parquet:
    input:
        rules.run_skani.output
    output:
        "Binning/{binner}/genome_similarities.parquet"
    run:
        import pandas as pd
        df = pd.read_table(input[0])

rule skani:
    input:
        "Binning/{binner}/genome_similarities.parquet".format(binner=config["final_binner"]),

