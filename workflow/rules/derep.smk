


rule skani:
    input:
        paths="Binning/{binner}/filtered_bins_paths.txt"
    output:
        "Intermediate/dereplication/{binner}_distance_matrix.txt"
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
        "skani sketch "
        " -l {input.paths} "
        " -o {output} "
        " -t {threads} "
        " --trace "
        " --sparse --ci "
        " --min-af 15"
        " &> {log} "



