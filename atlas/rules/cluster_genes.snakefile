
# this snakefile is to be executed in a single job of the master workflow in a cd-hit environment.
# It allows to eficiently sub-cluster genes whose product is part of the same protein cluster.
# The following parameters are expected in the config
#
# unclustered_dir
# clustered_dir
# coverage
# identity
import os


unclustered_dir= config['unclustered_dir']
clustered_dir = config['clustered_dir']



proteinIds, = glob_wildcards(os.path.join(unclustered_dir,"{proteinID}.fna"))


rule all:
    input:
        expand(os.path.join(clustered_dir,"{proteinID}.fna"), proteinID=proteinIds)


rule cluster_catalog:
    input:
        os.path.join(unclustered_dir,"{proteinID}.fna")
    output:
        os.path.join(clustered_dir,"{proteinID}.fna"),
        os.path.join(clustered_dir,"{proteinID}.clstr")
    params:
        prefix= lambda wc,output: os.path.splitext(output[0])[0],
    resources:
        mem=100 # can be reduced by comand line
    shell:
        """
            cd-hit-est -i {input} -T {threads} \
            -M {resources.mem}000 -o {params.prefix} \
            -c {config[identity]} -n 9  -d 0 \
            -aS {config[coverage]} -aL {config[coverage]}

            mv {params.prefix} {output[0]}
        """
