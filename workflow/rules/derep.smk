


rule run_skani:
    input:
        paths="Binning/{binner}/filtered_bins_paths.txt"
    output:
        temp("Intermediate/dereplication/{binner}_distance_matrix.txt")
    log:
        "logs/binning/{binner}/dereplication/skani_calculation.log"
    resources:
        mem_mb= config["mem"]*1000 ,
        time_min=60 * config["runtime"]["default"],
    params:
        #preset= "medium", # fast, medium or slow
        min_af=config["genome_dereplication"]["overlap"]*100,
        extra=""
    threads:
        config["threads"],
    conda:
        "../envs/skani.yaml"
    shell:
        "skani triangle "
        " {params.extra} "
        " -l {input.paths} "
        " -o {output} "
        " -t {threads} "
        " --trace "
        " --robust "
        " --sparse --ci "
        " --min-af {params.min_af} "
        " &> {log} "



rule skani_2_parquet:
    input:
        rules.run_skani.output
    output:
        "Binning/{binner}/genome_similarities.parquet"
    resources:
        mem_mb= config["mem"]*1000 ,
        time_min=60 * config["runtime"]["simplejob"],
    log:
        "logs/binning/{binner}/dereplication/skani_2_parquet.log"
    threads: 
        1
    run:

        try:
            import pandas as pd

            from utils.genome_dists import skani_column_dtypes


            import pandas as pd
            df = pd.read_table(input[0])

            from utils.io import simplify_path

            df = pd.read_table(input[0], usecols=list(skani_column_dtypes.keys()), dtype=skani_column_dtypes)

            df["Ref"] = df.Ref_file.cat.rename_categories(simplify_path)
            df["Query"] = df.Query_file.cat.rename_categories(simplify_path)

            df.to_parquet(output[0])
        
        except Exception as e:
            import traceback

            with open(log[0], "w") as logfile:
                traceback.print_exc(file=logfile)

            raise e
        

rule cluster_species:
    input:
        dist= "Binning/{binner}/genome_similarities.parquet",
        bin_info = "Binning/{binner}/filtered_bin_info.tsv"
    params:
        linkage_method = "average",
        pre_cluster_threshold = 0.925,
        threshold = config["genome_dereplication"]["ANI"]
    conda:
        "../envs/species_clustering.yaml"
    log:
        "logs/binning/{binner}/dereplication/species_clustering.log"
    output:
        cluster_file = "Binning/{binner}/bins2species.tsv"
    script:
        "../scripts/cluster_species.py"


rule skani:
    input:
        "Binning/{binner}/genome_similarities.parquet".format(binner=config["final_binner"]),
        "Binning/{binner}/bins2species.tsv".format(binner=config["final_binner"])

