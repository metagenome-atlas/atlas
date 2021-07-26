def parse_cd_hit_file(clstr_file):
    """

    >Cluster 0
    0   342nt, >S1_83_1... *
    1   342nt, >S2_82_1... at +/100.00%
    >Cluster 1
    0   339nt, >S1_61_1... *
    1   339nt, >S2_59_1... at +/100.00%


    """
    import numpy as np

    def parse_line(line):
        _, length, name, identity = (
            line.strip().replace("...", "\t").replace(", ", "\t").split("\t")
        )

        length = int(length.replace("nt", ""))
        name = name[1:]
        if "*" in identity:
            identity = np.nan
        else:
            identity = float(identity[identity.rfind("/") + 1 : identity.rfind("%")])

        return name, length, identity

    Clusters = []
    with open(clstr_file) as f:
        for line in f:
            if line[0] == ">":  # new cluster
                cluster = dict(elements=[], representative=None)
                Clusters.append(cluster)
            else:
                name, length, identity = parse_line(line)
                cluster["elements"].append((name, length, identity))
                if np.isnan(identity):
                    cluster["representative"] = name
    return Clusters


def write_cd_hit_clusters(Clusters, file_handle):
    for cluster in Clusters:
        for element in cluster["elements"]:
            file_handle.write(
                f"{element[0]}\t{element[1]}\t{element[2]}\t{cluster['representative']}\n"
            )


localrules:
    parse_clstr_files,
    rename_gene_clusters,


rule cluster_genes:
    input:
        fna_dir="Genecatalog/all_genes/predicted_genes.fna",
    output:
        temp("Genecatalog/representatives_of_clusters.fasta"),
        temp("Genecatalog/gene_catalog_oldnames.clstr"),
    conda:
        "%s/cd-hit.yaml" % CONDAENV
    log:
        "logs/Genecatalog/cluster_genes.log",
    threads: config.get("threads", 1)
    resources:
        mem=config["mem"],
    params:
        coverage=config["genecatalog"]["coverage"],
        identity=config["genecatalog"]["minid"],
        extra=config["genecatalog"]["extra"],
        prefix=lambda wc, output: os.path.splitext(output[1])[0],
    shell:
        """
        cd-hit-est -i {input} -T {threads} \
        -M {resources.mem}000 -o {params.prefix} \
        -c {params.identity} -n 9  -d 0 {params.extra} \
        -aS {params.coverage} -aL {params.coverage} &> {log}

        mv {params.prefix} {output[0]} 2>> {log}
        """


rule parse_clstr_files:
    input:
        clustered_dir="Genecatalog/gene_catalog_oldnames.clstr",
    output:
        temp("Genecatalog/clustering/orf2gene_oldnames.tsv"),
    run:
        with open(output[0], "w") as fout:
            fout.write(f"ORF\tLength\tIdentity\tGene\n")
            Clusters = parse_cd_hit_file(input[0])
            write_cd_hit_clusters(Clusters, fout)


rule rename_gene_clusters:
    input:
        orf2gene="Genecatalog/clustering/orf2gene_oldnames.tsv",
    output:
        orf2gene="Genecatalog/clustering/orf2gene.tsv.gz",
    run:
        import pandas as pd
        from Bio import SeqIO

        orf2gene = pd.read_csv(input.orf2gene, index_col=0, sep="\t")

        # rename gene repr to Gene0000XX

        gene_clusters_old_names = orf2gene["Gene"].unique()

        map_names = dict(
            zip(
                gene_clusters_old_names,
                utils.gen_names_for_range(len(gene_clusters_old_names), "Gene"),
            )
        )

        orf2gene["Gene"] = orf2gene["Gene"].map(map_names)
        orf2gene.to_csv(output.orf2gene, sep="\t", header=True, compression="gzip")
