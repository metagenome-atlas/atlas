

bin_quality_input_folder = "{sample}/binning/{binner}/bins"

## GUNC ###


#rule run_gunc:




##### BUSCO  #########

# rule run_busco:



##### checkM  #########
rule run_checkm_lineage_wf:
    input:
        touched_output="logs/checkm_init.txt",
        bins=bin_quality_input_folder,  # actualy path to fastas
    output:
        "{sample}/binning/{binner}/checkm/completeness.tsv",
        "{sample}/binning/{binner}/checkm/storage/tree/concatenated.fasta",
    params:
        output_dir=lambda wc, output: os.path.dirname(output[0]),
        tmpdir=config["tmpdir"],
    conda:
        "%s/checkm.yaml" % CONDAENV
    threads: config.get("threads", 1)
    log:
        "{sample}/logs/binning/{binner}/checkm.log",
    resources:
        time=config["runtime"]["long"],
        mem=config["large_mem"],
    benchmark:
        "logs/benchmarks/checkm_lineage_wf/{sample}_{binner}.tsv"
    shell:
        """
        rm -r {params.output_dir}
        checkm lineage_wf \
            --file {params.output_dir}/completeness.tsv \
            --tmpdir {params.tmpdir} \
            --tab_table \
            --quiet \
            --extension fasta \
            --threads {threads} \
            {input.bins} \
            {params.output_dir} &> {log}
        """


rule run_checkm_tree_qa:
    input:
        tree="{checkmfolder}/completeness.tsv",
    output:
        summary="{checkmfolder}/taxonomy.tsv",
    params:
        tree_dir=lambda wc, input: os.path.dirname(input.tree),
    conda:
        "%s/checkm.yaml" % CONDAENV
    threads: 1
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],
    shell:
        """
        checkm tree_qa \
           {params.tree_dir} \
           --out_format 2 \
           --file {output.summary}\
           --tab_table

        """


rule checkm_tetra:
    input:
        contigs=BINNING_CONTIGS,
    output:
        "{sample}/binning/{binner}/checkm/tetranucleotides.txt",
    log:
        "{sample}/logs/binning/{binner}/checkm/tetra.txt",
    conda:
        "%s/checkm.yaml" % CONDAENV
    threads: config.get("threads", 8)
    shell:
        """
        checkm tetra \
        --threads {threads} \
        {input.contigs} {output} 2> {log}
        """


rule checkm_outliers:
    input:
        tetra="{sample}/binning/{binner}/checkm/tetranucleotides.txt",
        bin_folder="{sample}/binning/{binner}/bins",
        checkm="{sample}/binning/{binner}/checkm/completeness.tsv",
    params:
        checkm_folder=lambda wc, input: os.path.dirname(input.checkm),
        report_type="any",
        treshold=95,  #reference distribution used to identify outliers; integer between 0 and 100 (default: 95)
    output:
        "{sample}/binning/{binner}/checkm/outliers.txt",
    log:
        "{sample}/logs/binning/{binner}/checkm/outliers.txt",
    conda:
        "%s/checkm.yaml" % CONDAENV
    threads: config.get("threads", 8)
    shell:
        """
        checkm outliers \
        --extension fasta \
        --distributions {params.treshold} \
        --report_type {params.report_type} \
        {params.checkm_folder} \
        {input.bin_folder} \
        {input.tetra} \
        {output} 2> {log}
        """
