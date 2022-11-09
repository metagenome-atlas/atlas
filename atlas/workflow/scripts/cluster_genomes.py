import pandas as pd

import scipy.spatial as sp
import scipy.cluster.hierarchy as hc
import numpy as np
from common import genome_pdist as gd
import networkx as nx



def threshold_based_clustering(Dist, threshold, linkage_method="average"):
    assert (threshold > 0.9) & (
        threshold < 1
    ), "threshold should be between 0.9 and 1 or 'auto', threshold was {threshold}"

    linkage = hc.linkage(sp.distance.squareform(Dist), method=linkage_method)

    labels = hc.fcluster(linkage, (1 - threshold), criterion="distance")


    return labels


if __name__ == "__main__":

    linkage_method = snakemake.params.linkage_method
    threshold = snakemake.params.threshold
    quality_score_formula = snakemake.config["quality_score"]

    Q = pd.read_csv(snakemake.input.genome_info,sep='\t',index_col=0)
    quality_score = Q.eval(quality_score_formula)

    assert (
        not quality_score.isnull().any()
    ), "I have NA quality values for thq quality score, it seems not all of the values defined in the quality_score_formula are presentfor all entries in tables/Genome_quality.tsv "



    # load genome distances
    M= gd.load_parquet(snakemake.input.dists)

    # genome distance to graph if ANI > 0.9
    G= gd.to_graph(M.query(f"Identity>=@threshold"))
    if hasattr(G,'selfloop_edges'):
        G.remove_edges_from(G.selfloop_edges())


    # prepare table for species number
    mag2Species = pd.DataFrame(index=Q.index, columns=["SpeciesNr", "Species"])
    mag2Species.index.name = "genome"

    last_species_nr=0

    for cc in nx.connected_components(G):


        Mcc = M.loc[( M.index.levels[0].intersection(cc),
                     M.index.levels[1].intersection(cc)
                     ),
                ]

        Dist = 1 - gd.pairewise2matrix(Mcc, fillna=Mcc.Identity.min())

        if threshold == "auto":
            Scores, labels = automatic_cluster_species(Dist, linkage_method=linkage_method)
        else:
            Scores, labels = threshold_based_clustering(
                Dist, threshold, linkage_method=linkage_method
            )

        # enter values of labels to species table
        mag2Species.loc[Dist.index, "SpeciesNr"] = labels +last_species_nr
        last_species_nr = mag2Species.SpeciesNr.max()


    missing_species = mag2Species.SpeciesNr.isnull()
    N_missing_species = sum(missing_species)
    mag2Species.loc[missing_species, "SpeciesNr"] = (
        np.arange(last_species_nr, last_species_nr + N_missing_species) + 1
    )

    Scores["N_clusters"] += N_missing_species
    Scores.to_csv(snakemake.output.scores, sep="\t", index=False)

    # assert (
    #     Scores.loc[Scores.Silhouette_score.idxmax(), "N_clusters"]
    #     == mag2Species.SpeciesNr.max()
    # ), "error in calculation of N species"
    print(f"Identified { mag2Species.SpeciesNr.max()} species")

    # create propper species names
    n_leading_zeros = len(str(max(labels)))
    format_int = "sp{:0" + str(n_leading_zeros) + "d}"
    mag2Species["Species"] = mag2Species.SpeciesNr.apply(format_int.format)

    # select representative
    mag2Species["Representative_Species"] = gd.best_genome_from_table(
        mag2Species.Species, quality_score
    )

    mag2Species.to_csv(snakemake.output.cluster_file, sep="\t")

    Q= Q.join(mag2Species)
    Q.to_csv(snakemake.input.genome_info,sep='\t')
