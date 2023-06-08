import os, sys
import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.DEBUG,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logging.captureWarnings(True)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logging.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


# Install exception handler
sys.excepthook = handle_exception

#### Begining of scripts


import pandas as pd

import scipy.spatial as sp
import scipy.cluster.hierarchy as hc
import numpy as np
from utils import genome_dist as gd
import networkx as nx


def get_float(value):
    "Enshure that value is [0-1]"

    assert value >= 0

    if value > 1:
        assert value <= 100, "it should be a percentage"
        logging.debug(f"Value {value} is > 1, I divede it with 100 to get a float")

        return value / 100
    else:
        return value


linkage_method = snakemake.params.linkage_method
pre_cluster_threshold = get_float(snakemake.params.pre_cluster_threshold)
threshold = get_float(snakemake.params.threshold)
min_aligned_fraction = get_float(snakemake.config["genome_dereplication"]["overlap"])

# verify ranges
gd.verify_expected_range(pre_cluster_threshold, 0.8, 1, "pre_cluster_threshold")
gd.verify_expected_range(threshold, 0.8, 1, "ANI cluster threshold")
gd.verify_expected_range(min_aligned_fraction, 0.1, 0.95, "min_aligned_fraction")


# load quality
Q = pd.read_csv(snakemake.input.bin_info, sep="\t", index_col=0)
# calculate quality score
logging.info("use quality score defined in bin info table")
quality_score = Q.quality_score

assert (
    not quality_score.isnull().any()
), "I have NA quality values for thq quality score, it seems not all of the values defined in the quality_score_formula are presentfor all entries in tables/Genome_quality.tsv "

logging.info("Load distances")
M = gd.load_skani(snakemake.input.dist)

# genome distance to graph
pre_clustering_criteria = (
    f"ANI >= {pre_cluster_threshold} & Align_fraction > {min_aligned_fraction}"
)

logging.info(f"Pre-cluster genomes with the if '{pre_clustering_criteria}'")
G = gd.to_graph(M.query(pre_clustering_criteria))

if hasattr(G, "selfloop_edges"):
    G.remove_edges_from(G.selfloop_edges())


# prepare table for species number
mag2Species = pd.DataFrame(index=Q.index, columns=["SpeciesNr", "Species"])
mag2Species.index.name = "genome"

last_species_nr = 0

n_pre_clusters = nx.connected.number_connected_components(G)
logging.info(f"Found {n_pre_clusters} pre-clusters, itterate over them.")
logging.debug(f"Cluster with threshold {threshold} and {linkage_method}-linkage method")
for i, cc in enumerate(nx.connected_components(G)):
    logging.info(f"Precluster {i+1}/{n_pre_clusters} with {len(cc)} genomes")

    Mcc = M.loc[
        (M.index.levels[0].intersection(cc), M.index.levels[1].intersection(cc)),
    ]

    labels = gd.group_species_linkage(
        Mcc.ANI, threshold=threshold, linkage_method=linkage_method
    )

    logging.debug(f"Got {labels.max()} species cluster for this pre-cluster.")

    # enter values of labels to species table
    mag2Species.loc[labels.index, "SpeciesNr"] = labels + last_species_nr
    last_species_nr = mag2Species.SpeciesNr.max()


missing_species = mag2Species.SpeciesNr.isnull()
N_missing_species = sum(missing_species)

logging.info(
    f"{N_missing_species} genomes were not part of a pre-cluster and are singletons."
)

mag2Species.loc[missing_species, "SpeciesNr"] = (
    np.arange(last_species_nr, last_species_nr + N_missing_species) + 1
)

logging.info(f"Identified { mag2Species.SpeciesNr.max()} species in total")

# create propper species names
n_leading_zeros = len(str(mag2Species.SpeciesNr.max()))
format_int = "sp{:0" + str(n_leading_zeros) + "d}"
mag2Species["Species"] = mag2Species.SpeciesNr.apply(format_int.format)


# select representative
logging.info("Select representative")
mag2Species["Representative_Species"] = gd.best_genome_from_table(
    mag2Species.Species, quality_score
)

mag2Species.to_csv(snakemake.output.cluster_file, sep="\t")

# Q= Q.join(mag2Species)
# Q.to_csv(snakemake.input.genome_info,sep='\t')
