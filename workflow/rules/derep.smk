


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
    threads: 
        1
    run:
        import pandas as pd

        skani_column_dtypes= {"Ref_file": "category",
                              "Query_file": "category",
                              "ANI": "float16",
                              "Align_fraction_ref": "float16",
                              "Align_fraction_query": "float16",
                              "ANI_5_percentile": "float16",
                              "ANI_95_percentile": "float16",
                              }# Ref_name        Query_name


        import pandas as pd
        from pathlib import Path
        df = pd.read_table(input[0])

        def simplify_filepath(path):
            "assumes single index are path of files, removes extesnion and dirname"

            path = Path(path)
            
            assert len(path.suffixes) ==1,f"{path} has multiple extensions. gzip is not yet supported"
            return path.stem

        df = pd.read_table(input[0], use_cols=list(skani_column_dtypes.keys()), dtypes=skani_column_dtypes)

        df["Ref"] = df.Ref_file.cat.rename_categories(simplify_filepath)
        df["Query"] = df.Query_file.cat.rename_categories(simplify_filepath)

        df.to_parquet(output[0])
        

rule skani:
    input:
        "Binning/{binner}/genome_similarities.parquet".format(binner=config["final_binner"]),

