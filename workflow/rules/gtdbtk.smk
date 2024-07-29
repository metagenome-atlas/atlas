gtdb_dir = "genomes/taxonomy/gtdb"


rule identify:
    input:
        flag=rules.extract_gtdb.output,
        genes_flag="genomes/annotations/genes/predicted",
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
        extension="faa",
        gene_dir=lambda wc, input: os.path.abspath(os.path.dirname(input.genes_flag)),
    shell:
        'export GTDBTK_DATA_PATH="{GTDBTK_DATA_PATH}" ; '
        "gtdbtk identify "
        "--genes --genome_dir {params.gene_dir} "
        " --out_dir {params.outdir} "
        "--extension {params.extension} "
        "--cpus {threads} &> {log[0]}"


checkpoint align:
    input:
        f"{gtdb_dir}/identify",
    output:
        directory(f"{gtdb_dir}/align"),
    threads: config["threads"]
    resources:
        mem_mb=config["large_mem"] * 1000,
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk/align.txt",
        f"{gtdb_dir}/gtdbtk.log",
    params:
        outdir=gtdb_dir,
    shell:
        'export GTDBTK_DATA_PATH="{GTDBTK_DATA_PATH}" ; '
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
        mem_mb=config["large_mem"] * 1000,
        time_min=60 * config["runtime"]["long"],
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/taxonomy/gtdbtk/classify.txt",
        f"{gtdb_dir}/gtdbtk.log",
    params:
        outdir=gtdb_dir,
        extension="fasta",
        mashdir=Path(GTDBTK_DATA_PATH) / "mash_db",
    shell:
        'export GTDBTK_DATA_PATH="{GTDBTK_DATA_PATH}" ; '
        "gtdbtk classify --genome_dir {input.genome_dir} --align_dir {params.outdir} "
        " --mash_db {params.mashdir} "
        "--out_dir {params.outdir} "
        " --tmpdir {resources.tmpdir} "
        "--extension {params.extension} "
        "--cpus {threads} &> {log[0]}"


rule combine_taxonomy:
    input:
        folder=f"{gtdb_dir}/classify",
    output:
        combined=f"{gtdb_dir}/gtdbtk.combined.summary.tsv",
        taxonomy="genomes/taxonomy/gtdb_taxonomy.tsv",
    log:
        "logs/taxonomy/gtdbtk/combine.txt",
    script:
        "../scripts/combine_taxonomy.py"


rule build_tree:
    input:
        f"{gtdb_dir}/align/{{msa}}.user_msa.fasta.gz",
    output:
        temp("genomes/taxonomy/gtdb/{msa}.unrooted.tree"),
    log:
        "logs/genomes/tree/{msa}.log",
        "logs/genomes/tree/{msa}.err",
    threads: max(config["threads"], 3)
    params:
        outdir=lambda wc, output: Path(output[0]).parent,
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        'export GTDBTK_DATA_PATH="{GTDBTK_DATA_PATH}" ; '
        "gtdbtk infer --msa_file {input} "
        " --out_dir {params.outdir} "
        " --prefix {wildcards.msa} "
        " --cpus {threads} "
        "--tmpdir {resources.tmpdir} > {log[0]} 2> {log[1]}"


localrules:
    root_tree,


rule root_tree:
    input:
        tree=rules.build_tree.output[0],
    wildcard_constraints:
        msa="((?!unrooted).)*",
    output:
        tree="genomes/tree/{msa}.nwk",
    conda:
        "../envs/tree.yaml"
    threads: 1
    resources:
        mem_mb=config["simplejob_mem"] * 1000,
        ttime_min=60 * config["runtime"]["simplejob"],
    log:
        "logs/genomes/tree/root_tree_{msa}.log",
    script:
        "../scripts/root_tree.py"


def all_gtdb_trees_input(wildcards):
    dir = checkpoints.align.get().output[0]

    domains = glob_wildcards(f"{dir}/gtdbtk.{{domain}}.user_msa.fasta.gz").domain

    return expand("genomes/tree/gtdbtk.{domain}.nwk", domain=domains)


rule all_gtdb_trees:
    input:
        all_gtdb_trees_input,
    output:
        touch("genomes/tree/finished_gtdb_trees"),
