


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
    #params:
        #preset= "medium" # fast, medium or slow
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

            skani_column_dtypes= {"Ref_file": "category",
                                "Query_file": "category",
                                "ANI": float,
                                "Align_fraction_ref": float,
                                "Align_fraction_query": float,
                                "ANI_5_percentile": float,
                                "ANI_95_percentile": float,
                                }# Ref_name        Query_name


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
        

rule skani:
    input:
        "Binning/{binner}/genome_similarities.parquet".format(binner=config["final_binner"]),

