DBDIR = config["database_dir"]


def get_dram_config(wildcards):
    return config.get("dram_config_file", f"{DBDIR}/DRAM.config")


rule dram_download:
    output:
        dbdir=directory(f"{DBDIR}/Dram/"),
        config=f"{DBDIR}/DRAM.config",
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["default"],
    log:
        "log/dram/download_dram.log",
    benchmark:
        "log/benchmarks/dram/download_dram.tsv"
    conda:
        "../envs/dram.yaml"
    shell:
        " DRAM-setup.py prepare_databases "
        " --output_dir {output.dbdir} "
        " --threads {threads} "
        " --verbose "
        " --skip_uniref "
        " &> {log} "
        " ; "
        " DRAM-setup.py export_config --output_file {output.config}"


localrules:
    DRAM_set_db_loc,
    get_lists_of_genomes,


rule DRAM_set_db_loc:
    input:
        get_dram_config,
    output:
        touch("logs/dram_config_imported"),
    threads: 1
    conda:
        "../envs/dram.yaml"
    shell:
        "DRAM-setup.py import_config --config_loc {input}"


checkpoint get_lists_of_genomes:
    input:
        get_genome_folder,
    output:
        directory(temp("annotations/dram/genome_lists")),
    run:
        from glob import glob

        all_fasta_files = glob(input[0] + "/*.fasta")

        if len(all_fasta_files) == 0:
            raise Exception(f"No genome fasta files found in folder {input[0]}")

        os.makedirs(output[0])

        N = 10  # N subsets
        for i in range(0, len(all_fasta_files) // N + 1):

            with open(output[0] + f"/subset_{i+1}.txt", "w") as outf:
                outf.write(" ".join(all_fasta_files[i * N : (i + 1) * N]))


def dram_annotate_input(wildcards, input):

    fasta_files = open(input.genome_list).read().strip().split()

    assert len(fasta_files) > 0

    return " ".join(f"--input_fasta {fasta}" for fasta in fasta_files)


rule DRAM_annotate:
    input:
        genome_folder=get_genome_folder,
        genome_list="annotations/dram/genome_lists/{subset}.txt",
        flag=rules.DRAM_set_db_loc.output,
    output:
        outdir=directory("annotations/dram/intermediate_files/{subset}"),
    threads: config["threads"]
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["default"],
    conda:
        "../envs/dram.yaml"
    shadow:
        "shallow"
    params:
        gtdb_file="gtdbtk.bac120.summary.tsv",
        input=dram_annotate_input,
    log:
        "log/dram/run_dram/{subset}.log",
    benchmark:
        "log/benchmarks/dram/run_dram/{subset}.tsv"
    shell:
        " DRAM.py annotate "
        "  {params.input} "
        " --output_dir {output.outdir} "
        " --prodigal_mode single "
        " --threads {threads} "
        " --verbose &> {log}"


def get_all_dram(wildcards):

    genome_lists_folder = checkpoints.get_lists_of_genomes.get().output[0]

    all_subsets = glob_wildcards(
        os.path.join(genome_lists_folder, "{subset}.txt")
    ).subset

    if len(all_subsets) == 0:

        raise Exception(
            f"No genome list file is found in {genome_lists_folder}\n"
            "Remove this empty folder and start anew."
        )

    return expand(rules.DRAM_annotate.output.outdir, subset=all_subsets)


DRAM_ANNOTATON_FILES = ["annotations.tsv", "rrnas.tsv", "trnas.tsv"]


localrules:
    concat_annotations,


rule concat_annotations:
    input:
        get_all_dram,
    output:
        expand("annotations/dram/{annotation}", annotation=DRAM_ANNOTATON_FILES),
    resources:
        time=config["runtime"]["default"],
    #     mem = config['mem']
    run:
        # from utils import io
        for i, annotation_file in enumerate(DRAM_ANNOTATON_FILES):

            input_files = [
                os.path.join(dram_folder, annotation_file) for dram_folder in input
            ]

            # drop files that don't exist for rrna and trna
            if not i == 0:
                input_files = [f for f in input_files if os.path.exists(f)]

            shell(f"head -n1 {input_files[0]} > {output[i]} ")
            for f in input_files:
                shell(f"tail -n+2 {f} >> {output[i]}")

        # io.pandas_concat(input_files, output[i],sep='\t',index_col=0, axis=0)



rule DRAM_destill:
    input:
        expand("annotations/dram/{annotation}", annotation=DRAM_ANNOTATON_FILES),
        flag=rules.DRAM_set_db_loc.output,
    output:
        outdir=directory("annotations/dram/distil"),
    threads: 1
    resources:
        mem=config['simplejob_mem'],
        time=config["simplejob_time"]
    conda:
        "../envs/dram.yaml"
    log:
        "log/dram/distil.log",
    shell:
        " DRAM.py distill "
        " --input_file {input[0]}"
        " --rrna_path {input[1]}"
        " --trna_path {input[2]}"
        " --output_dir {output} "
        "  &> {log}"


rule get_all_modules:
    input:
        "annotations/dram/annotations.tsv",
    output:
        "annotations/dram/kegg_modules.tsv",
    threads: 1
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["default"],
    conda:
        "../envs/dram.yaml"
    log:
        "log/dram/get_all_modules.log",
    script:
        "../scripts/DRAM_get_all_modules.py"


rule dram:
    input:
        "annotations/dram/distil",
        "annotations/dram/kegg_modules.tsv",
