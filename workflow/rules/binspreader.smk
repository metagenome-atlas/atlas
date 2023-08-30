
rule rename_paths:
    input:
        paths="{{sample}}/assembly/{sequences}.paths".format(
            sequences="scaffolds" if config["spades_use_scaffolds"] else "contigs"
            ),
        filtered_contigs = get_assembly
    output:
        "{sample}/assembly/renamed_assembly.paths"
    run:

        with open(output[0],"w") as out_handle, open(input.paths) as in_handle:
            for line in in_handle:
                if line.startswith("NODE"):

                    # expect"NODE_1_length_196862_cov_1048.163282"
                    _,nr,_ = line.split("_",2)

                    out_handle.write(f"{wildcards.sample}_{nr}\n")
                else:
                    out_handle.write(line)





rule prepare_yaml_for_binspreader:
    input:
        reads=get_quality_controlled_reads
    output:
        temp("Intermediate/binspreader/dataset_yamls/{sample}.yaml")
    run:
        assert len(input.reads)>=2

        import yaml
        from os.path import abspath

        data = [
            {
                'orientation': 'fr',
                'type': 'paired-end',
                'right reads': [abspath(input.reads[0])],
                'left reads': [abspath(input.reads[1])]
            }
        ]

        with open(output[0], 'w') as file:
            yaml.dump(data, file, default_style='"', default_flow_style=False)

            


rule binspreader:
    input:
        #get_quality_controlled_reads,
        graph="{sample}/assembly/assembly_graph_after_simplification.gfa",
        #yaml = ancient(rules.prepare_yaml_for_binspreader.output[0]),
        paths = "{sample}/assembly/renamed_assembly.paths",
        bins = "{sample}/binning/{binner_before}/cluster_attribution.tsv",
    output:
        directory("Intermediate/binspreader/after_{binner_before}/outputs/{sample}")
    log:
        "logs/binspreader/after_{binner_before}/run_binspreader/{sample}.log"
    conda:
        "../envs/binspreader.yaml"
    threads:
        config["threads"]
    shell:
        "bin-refine "
        " {input.graph} " 
        " {input.bins} "
        " {output} " 
        " --paths {input.paths} "
        " --no-unbinned-bin "
        "  -Rcorr "
        " -m " # allow multiple bin assignment
        " --tmp-dir {resources.tmpdir}/binspreader_{wildcards.sample}_{wildcards.binner_before} "
        " &> {log} "

        #"  " # only part of graph is binned
        #" --reads "
        #" --dataset {input.yaml} "
        # " --sparse-propagation"


rule parse_binspreader:
    input:
        rules.binspreader.output
    output:
        "{sample}/binning/binsp{binner_before}/cluster_attribution.tsv"
    run:

        input_file = Path(input[0])/"binning.tsv"

        with open(output[0],"w") as out_handle, open(input_file) as in_handle:
            for line in in_handle:
                # contig tab binname1 tab binname2 ...

                contig,bin_names = line.strip().split('\t',1)

                for bin in bin_names.split():
                    new_bin_name = bin.replace(wildcards.binner_before,"binsp"+wildcards.binner_before)

                    out_handle.write(f"{contig}\t{new_bin_name}\n")