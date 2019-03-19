


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

# rule root_tree:
#     input:
#         tree"genomes/checkm/storage/tree/concatenated.fasta",
#         taxonomy="genomes/taxonomy/taxonomy.tsv"
#     output:
#         "genomes/tree.nwk"
#     conda:
#         "%s/tree.yaml" % CONDAENV
#     threads:
#         1
#     run:
