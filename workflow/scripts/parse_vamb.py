import os, sys
import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


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


import pandas as pd
from pathlib import Path
from utils.utils import gen_names_for_range
from utils.fasta import parse_fasta_headers


fasta_extension = snakemake.params.fasta_extension
separator = snakemake.params.separator

# path with {sample} to replace
cluster_output_path = snakemake.params.output_path


# cluster.tsv.gz file for all samples of all bingroups
output_culsters = snakemake.output.renamed_clusters


all_clusters = []


for i in range(len(list(snakemake.input))):
    vamb_folder = Path(snakemake.input[i])

    bingroup = snakemake.params.bingroups[i]

    logging.info(f"Parse vamb output for bingroup {bingroup}")

    # path to the bins folder
    bin_dir = vamb_folder / "bins"
    vamb_cluster_file = vamb_folder / "clusters.tsv"

    # Get a list of binds that are big enough. Not all bins in the vamb_cluster_file, pass the size filter
    big_bins = []

    for file in os.listdir(bin_dir):
        bin_name, extension = os.path.splitext(file)

        logging.debug(f"Found file {bin_name} with extension {extension}")

        if extension == fasta_extension:
            big_bins.append(bin_name)

    logging.info(
        f"Found {len(big_bins)} bins created by Vamb (above size limit)\n"
        f"E.g. {big_bins[:5]}"
    )

    logging.info(f"Load vamb cluster file {vamb_cluster_file}")
    clusters_contigs = pd.read_table(vamb_cluster_file, header=None)
    clusters_contigs.columns = ["OriginalName", "Contig"]

    # split contigs by separator. This is mainly done for compatibility with SemiBin
    clusters = clusters_contigs.Contig.str.rsplit(separator, n=1, expand=True)
    clusters.columns = ["Sample", "Contig"]

    # get number of BinID given by vamb, prefix with bingroup
    clusters["BinId"] = (
        bingroup
        + clusters_contigs.OriginalName.str.rsplit(separator, n=1, expand=True)[1]
    )

    # Add information if the bin is large enough
    clusters["OriginalName"] = clusters_contigs.OriginalName
    clusters["Large_enough"] = clusters.OriginalName.isin(big_bins)

    # Add information about the bingroup
    clusters["BinGroup"] = bingroup

    all_clusters.append(clusters)

    del clusters_contigs


logging.info(f"Concatenate clusters of all bingroups")
clusters = pd.concat(all_clusters, axis=0)


n_bins = (
    clusters.query("Large_enough").groupby(["BinGroup", "Sample"])["BinId"].nunique()
)
logging.info(
    f"Number of bins per sample and bingroup passing the size filter:\n{n_bins}"
)


clusters["SampleBin"] = clusters.Sample + "_vamb_" + clusters.BinId
clusters.loc[~clusters.Large_enough, "SampleBin"] = ""


logging.info(f"Write reformated table to {output_culsters}")
clusters.to_csv(output_culsters, sep="\t", index=False)

# filter for following
clusters = clusters.query("Large_enough")

logging.info(f"Write cluster_attribution for samples")
for sample, cl in clusters.groupby("Sample"):
    sample_output_path = cluster_output_path.format(sample=sample)

    logging.debug(f"Write file {sample_output_path}")
    cl[["Contig", "SampleBin"]].to_csv(
        sample_output_path, sep="\t", index=False, header=False
    )


samples_without_bins = set(snakemake.params.samples).difference(set(clusters.Sample))

if len(samples_without_bins) > 0:
    logging.warning(
        "The following samples did't yield bins, I add longest contig to make the pipline continue:\n"
        + "\n".join(samples_without_bins)
    )

    for sample in samples_without_bins:
        sample_output_path = cluster_output_path.format(sample=sample)
        with open(sample_output_path, "w") as fout:
            fout.write(f"{sample}_1\t{sample}_vamb_1\n")
