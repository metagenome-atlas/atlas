


rule fasttree:
    input:
        "genomes/checkm/storage/tree/concatenated.fasta"
    output:
        "genomes/checkm/storage/tree/fasttree.nwk"
    log:
        "logs/genomes/tree/FastTree.log"
    threads:
        config['threads']
    conda:
        "%s/tree.yaml" % CONDAENV
    shell:
        "export OMP_NUM_THREADS={threads}; "
        "FastTree -log {log} {input} > {output} "

localrules: root_and_plot_tree
rule root_and_plot_tree:
    input:
        tree="genomes/checkm/storage/tree/fasttree.nwk",
        taxonomy="genomes/checkm/taxonomy.tsv"
    output:
        tree="genomes/tree/tree.nwk",
        svg="genomes/tree/tree.svg",
        pdf="genomes/tree/tree.pdf"
    conda:
        "%s/tree.yaml" % CONDAENV
    threads:
        1
    script:
        "../scripts/utils/tree.py"
