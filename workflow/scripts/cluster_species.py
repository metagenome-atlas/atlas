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

#### Beginning of scripts


import pandas as pd

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
Q.Additional_Notes = Q.Additional_Notes.fillna("").astype(str)


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
genomes_to_drop = []

last_species_nr = 1  # start at 1


n_pre_clusters = nx.connected.number_connected_components(G)
logging.info(f"Found {n_pre_clusters} pre-clusters, iterate over them.")
logging.debug(f"Cluster with threshold {threshold} and {linkage_method}-linkage method")
for i, cc in enumerate(nx.connected_components(G)):
    logging.info(f"Precluster {i+1}/{n_pre_clusters} with {len(cc)} genomes")

    Qcc = Q.loc[list(cc)]

    # check translation table
    freq = Qcc["Translation_Table_Used"].value_counts()

    if freq.shape[0] > 1:
        logging.info(
            "Not all genomes use the same translation table,"
            "drop genomes that don't use main translation table."
        )
        logging.info(freq)

        main_tranlation_table = freq.index[0]

        drop_genomes = Qcc.query(
            "Translation_Table_Used != @main_tranlation_table"
        ).index

        cc = cc - set(drop_genomes)
        Qcc = Qcc.loc[list(cc)]
        genomes_to_drop += list(drop_genomes)
        logging.info(f"Drop {len(drop_genomes) } genomes, keep ({len(cc)})")

    ## Check that the same completeness model is used for all

    freq = Qcc["Completeness_Model_Used"].value_counts()
    if freq.shape[0] > 1:
        logging.info(
            "Not all genomes use the same completeness model. Recalibrate completeness."
        )

        logging.info(freq)

        # genomes that don't use specific model
        non_specific = Qcc.index[
            ~Qcc.Completeness_Model_Used.str.contains("Specific Model")
        ]

        logging.debug(
            f"{len(non_specific)} genomes use general completeness model. Recalibrate completeness and quality score to use lower value"
        )

        logging.debug(
            Qcc.loc[
                non_specific,
                ["Completeness_General", "Completeness_Specific", "Contamination"],
            ]
        )

        Qcc.loc[non_specific, "Completeness"] = Qcc.loc[
            non_specific,
            [
                "Completeness_General",
                "Completeness_Specific",
            ],
        ].min(axis=1)

        # add note
        Q.loc[
            non_specific, "Additional_Notes"
        ] += "Completeness was re-calibrated based on Completeness model used in all genomes of the species."

        # transfer to main quality
        Q.loc[list(cc), "Completeness"] = Qcc.loc[list(cc), "Completeness"]

        # drop low quality genomes

        logging.info("Drop low quality genomes according to filter criteria")

        try:
            filter_criteria = snakemake.config["genome_filter_criteria"]
            drop_genomes = Qcc.index.difference(Qcc.query(filter_criteria).index)

        except Exception as e:
            logging.error("Cannot filter low quality genomes")
            logging.exception(e)

            drop_genomes = []

        if len(drop_genomes) > 0:
            cc = cc - set(drop_genomes)
            logging.info(
                f"Drop {len(drop_genomes) } with too low quality genomes, keep {len(cc)}"
            )

            Qcc = Qcc.loc[list(cc)]
            genomes_to_drop += list(drop_genomes)

    if len(cc) <= 1:
        logging.info(
            "I am left with {len(cc)} genomes in this pre-cluster. No need to cluster."
        )
    else:
        # subset dist matrix
        Mcc = M.loc[
            (M.index.levels[0].intersection(cc), M.index.levels[1].intersection(cc)),
        ]

        # Cluster species
        labels = gd.group_species_linkage(
            Mcc.ANI, threshold=threshold, linkage_method=linkage_method
        )

        logging.debug(f"Got {labels.max()} species cluster for this pre-cluster.")

        # enter values of labels to species table
        mag2Species.loc[labels.index, "SpeciesNr"] = labels + last_species_nr
        last_species_nr = mag2Species.SpeciesNr.max()


mag2Species.drop(genomes_to_drop, inplace=True)
Q.drop(genomes_to_drop, inplace=True)


missing_species = mag2Species.index[mag2Species.SpeciesNr.isnull()]
N_missing_species = len(missing_species)

logging.info(
    f"{N_missing_species} genomes were not part of a pre-cluster and are singleton-species."
)

Q.loc[missing_species, "Additional_Notes"] += " Singleton species"

mag2Species.loc[missing_species, "SpeciesNr"] = (
    np.arange(last_species_nr, last_species_nr + N_missing_species) + 1
)


n_species = mag2Species.SpeciesNr.unique().shape[0]
logging.info(f"Identified {n_species } species in total")

# create proper species names
n_leading_zeros = len(str(mag2Species.SpeciesNr.max()))
format_int = "sp{:0" + str(n_leading_zeros) + "d}"
mag2Species["Species"] = mag2Species.SpeciesNr.apply(format_int.format)


# calculate quality score


logging.info("Define Quality score defined as Completeness - 5x Contamination")
# recalculate quality score as some completeness might be recalibrated.
Q.eval("Quality_score = Completeness - 5* Contamination", inplace=True)
quality_score = Q.Quality_score

assert (
    not quality_score.isnull().any()
), "I have NA quality values for the quality score, it seems not all of the values defined in the quality_score_formula are presentfor all entries in tables/Genome_quality.tsv "


# select representative
logging.info("Select representative")
mag2Species["Representative"] = gd.best_genome_from_table(
    mag2Species.Species, quality_score
)

mag2Species.to_csv(snakemake.output.bins2species, sep="\t")

# mag2Species = mag2Species.join(Q)
Q.to_csv(snakemake.output.bin_info, sep="\t")
