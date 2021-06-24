gtdb_dir = "genomes/taxonomy/gtdb"


rule identify:
    input:
        dir=genome_dir,
        flag=rules.download_gtdb.output,
    output:
        directory(f"{gtdb_dir}/identify"),
    threads: config["threads"]
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk/identify.txt",
        f"{gtdb_dir}/gtdbtk.log",
    params:
        outdir=gtdb_dir,
        extension="fasta",
    shell:
        "GTDBTK_DATA_PATH={GTDBTK_DATA_PATH} ; "
        "gtdbtk identify --genome_dir {input.dir} --out_dir {params.outdir} "
        "--extension {params.extension} "
        "--cpus {threads} &> {log[0]}"


checkpoint align:
    input:
        f"{gtdb_dir}/identify",
    output:
        directory(f"{gtdb_dir}/align"),
    threads: config["threads"]
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk/align.txt",
        f"{gtdb_dir}/gtdbtk.log",
    params:
        outdir=gtdb_dir,
    shell:
        "GTDBTK_DATA_PATH={GTDBTK_DATA_PATH} ;  "
        "gtdbtk align --identify_dir {params.outdir} --out_dir {params.outdir} "
        "--cpus {threads} &> {log[0]}"


rule classify:
    input:
        rules.align.output,
        genome_dir=genome_dir,
    output:
        directory(f"{gtdb_dir}/classify"),
    threads: config["threads"]  #pplacer needs much memory for not many threads
    resources:
        mem=config["large_mem"],
        time=config["runtime"]["long"],
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk/classify.txt",
        f"{gtdb_dir}/gtdbtk.log",
    params:
        outdir=gtdb_dir,
        extension="fasta",
    shell:
        "GTDBTK_DATA_PATH={GTDBTK_DATA_PATH} ; "
        "gtdbtk classify --genome_dir {input.genome_dir} --align_dir {params.outdir} "
        "--out_dir {params.outdir} "
        "--extension {params.extension} "
        "--cpus {threads} &> {log[0]}"


msa_paths = {
    "checkm": "genomes/checkm/storage/tree/concatenated.fasta",
    "gtdbtk.bac120": f"{gtdb_dir}/align/gtdbtk.bac120.user_msa.fasta",
    "gtdbtk.ar122": f"{gtdb_dir}/align/gtdbtk.ar122.user_msa.fasta",
}


rule fasttree:
    input:
        lambda wildcards: msa_paths[wildcards.msa],
    output:
        temp("genomes/tree/{msa}.unrooted.nwk"),
    log:
        "logs/genomes/tree/FastTree_{msa}.log",
    threads: max(config["threads"], 3)
    conda:
        "%s/tree.yaml" % CONDAENV
    shell:
        "export OMP_NUM_THREADS={threads}; "
        "FastTree -log {log} {input} > {output} "


localrules:
    root_tree,


rule root_tree:
    input:
        tree="genomes/tree/{msa}.unrooted.nwk",
    wildcard_constraints:
        msa="((?!unrooted).)*",
    output:
        tree="genomes/tree/{msa}.nwk",
    conda:
        "%s/tree.yaml" % CONDAENV
    threads: 1
    log:
        "logs/genomes/tree/root_tree_{msa}.log",
    script:
        "../scripts/root_tree.py"


def all_gtdb_trees_input(wildcards):
    dir = checkpoints.align.get().output[0]

    domains = glob_wildcards(f"{dir}/gtdbtk.{{domain}}.user_msa.fasta").domain

    return expand("genomes/tree/gtdbtk.{domain}.nwk", domain=domains)


rule all_gtdb_trees:
    input:
        all_gtdb_trees_input,
    output:
        touch("genomes/tree/finished_gtdb_trees"),
