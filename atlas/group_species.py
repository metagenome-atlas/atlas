import pandas as pd

import scipy.spatial as sp
import scipy.cluster.hierarchy as hc
from sklearn.metrics import silhouette_score
import numpy as np
from common import genome_pdist as gd


def automatic_cluster_species(
    Dist, seed_tresholds=[0.92, 0.97], linkage_method="average"
):

    linkage = hc.linkage(sp.distance.squareform(Dist), method=linkage_method)

    def get_Nclusters(treshold):
        labels = hc.fcluster(linkage, (1 - treshold), criterion="distance")
        return max(labels)

    N_range = [get_Nclusters(t) for t in seed_tresholds]

    assert (N_range[1] - N_range[0]) < 60, "Need to evaluate more than 60 tresholds"

    assert ~np.isnan(N_range).any(), "N range is not defined"

    Scores = gd.evaluate_clusters_range(
        np.arange(min(N_range), max(N_range) + 1), Dist, linkage_method=linkage_method
    )

    if N_range[0] == N_range[1]:
        labels = hc.fcluster(linkage, (1 - seed_tresholds[0]), criterion="distance")
    else:

        N_species = Scores.Silhouette_score.idxmax()
        labels = hc.fcluster(linkage, N_species, criterion="maxclust")

    return Scores, labels


def treshold_based_clustering(Dist, treshold, linkage_method="average"):
    assert (treshold > 0.9) & (
        treshold < 1
    ), "treshold should be between 0.9 and 1 or 'auto', treshold was {treshold}"
    linkage = hc.linkage(sp.distance.squareform(Dist), method=linkage_method)
    labels = hc.fcluster(linkage, (1 - treshold), criterion="distance")
    Scores = gd.evaluate_clusters_tresholds(
        [treshold], Dist, linkage_method=linkage_method
    )
    return Scores, labels


if __name__ == "__main__":

    linkage_method = snakemake.params.linkage_method
    treshold = snakemake.params.treshold
    quality_score_formula = snakemake.config["quality_score"]

    Q = gd.load_quality(snakemake.input.quality)
    quality_score = Q.eval(quality_score_formula)

    M = gd.load_mummer(snakemake.input.dists)
    Dist = 1 - gd.pairewise2matrix(M, fillna=0.9)

    if treshold == "auto":
        Scores, labels = automatic_cluster_species(Dist, linkage_method=linkage_method)
    else:
        Scores, labels = treshold_based_clustering(
            Dist, treshold, linkage_method=linkage_method
        )

    Scores.to_csv(snakemake.output.scores, sep="\t")

    mag2Species = pd.DataFrame(index=Q.index, columns=["SpeciesNr", "Species"])
    mag2Species.index.name = "genome"
    mag2Species.loc[Dist.index, "SpeciesNr"] = labels

    speciesNr = labels.max()
    missing_species = mag2Species.SpeciesNr.isnull()
    mag2Species.loc[missing_species, "SpeciesNr"] = np.arange(
        speciesNr + 1, speciesNr + 1 + missing_species.sum()
    )

    print(f"Identified { mag2Species.SpeciesNr.max()} species")

    n_leading_zeros = len(str(max(labels)))
    format_int = "sp{:0" + str(n_leading_zeros) + "d}"
    mag2Species["Species"] = mag2Species.SpeciesNr.apply(format_int.format)

    mag2Species["Representative_Species"] = gd.best_genome_from_table(
        mag2Species.Species, quality_score
    )

    mag2Species.to_csv(snakemake.output.cluster_file, sep="\t")
