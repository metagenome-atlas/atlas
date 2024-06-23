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
        mem_mb=config["mem"] * 1000,
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
        temp("Genecatalog/orf2gene_oldnames.tsv"),
    run:
        with open(output[0], "w") as fout:
            fout.write(f"ORF\tLength\tIdentity\tRepresentative\n")
            Clusters = parse_cd_hit_file(input[0])
            write_cd_hit_clusters(Clusters, fout)


rule generate_orf_info:
    input:
        cluster_attribution="Genecatalog/orf2gene_oldnames.tsv",
    output:
        cluster_attribution="Genecatalog/clustering/orf_info.parquet",
        rep2genenr="Genecatalog/clustering/representative2genenr.tsv",
    threads: 1
    run:
        import pandas as pd
        import numpy as np

        from utils import gene_scripts

        # cd hit format ORF\tLength\tIdentity\tRepresentative\n
        orf2gene = pd.read_csv(input.orf2gene, sep="\t")

        # rename gene repr to Gene0000XX

        # split orf names in sample, contig_nr, and orf_nr
        orf_info = gene_scripts.split_orf_to_index(orf2gene.ORF)

        # rename representative

        representative_names = orf2gene.Representative.unique()

        map_names = pd.Series(
            index=representative_names,
            data=np.arange(1, len(representative_names) + 1, dtype=np.uint),
        )


        orf_info["GeneNr"] = orf2gene.Representative.map(map_names)


        orf_info.to_parquet(output.cluster_attribution)


        # Save name of representatives
        map_names.index.name = "Representative"
        map_names.name = "GeneNr"
        map_names.to_csv(output.rep2genenr, sep="\t")
