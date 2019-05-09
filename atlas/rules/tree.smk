


rule fasttree:
    input:
        "genomes/checkm/storage/tree/concatenated.fasta"
    output:
        "genomes/checkm/storage/tree/fasttree.nwk"
    log:
        "logs/genomes/tree/FastTree.log"
    threads:
        max(config['threads'],3)
    conda:
        "%s/tree.yaml" % CONDAENV
    shell:
        "export OMP_NUM_THREADS={threads}; "
        "FastTree -log {log} {input} > {output} "

localrules: root_tree
rule root_tree:
    input:
        tree="genomes/checkm/storage/tree/fasttree.nwk",
        taxonomy="genomes/checkm/taxonomy.tsv"
    output:
        tree="genomes/tree/tree.nwk",
    conda:
        "%s/tree.yaml" % CONDAENV
    threads:
        1
    script:
        "../scripts/utils/tree.py"
