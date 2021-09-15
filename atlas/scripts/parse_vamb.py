import os, sys
import logging, traceback

logging.basicConfig(
    #    filename=snakemake.log[0], # HACK
    #    level=logging.INFO,
    level=logging.DEBUG,
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
from utils.utils import gen_names_for_range
from utils.fasta import parse_fasta_headers

bin_dir = os.path.join(snakemake.input[0], "bins")
fasta_extension = snakemake.params.fasta_extension
separator = snakemake.params.separator

cluster_output_path = snakemake.params.output_path

vamb_cluster_file = os.path.join(snakemake.input[0], "clusters.tsv")
output_culsters = snakemake.output.renamed_clusters


big_bins = []

for file in os.listdir(bin_dir):
    binname, extension = os.path.split(file)
    if extension == fasta_extension:
        bins.append(big_bins)


logging.info(f"Load vamb cluster file {vamb_cluster_file}")
clusters = pd.read_table(vamb_cluster_file, header=None)

clusters_contigs.columns = ["OriginalName", "Contig"]


clusters = clusters_contigs.Contig.str.rsplit(separator, n=1, expand=True)
clusters.columns = ["Sample", "Contig"]

clusters["BinID"] = clusters_contigs.OriginalName.str.rsplit(
    separator, n=1, expand=True
)[1]
clusters["OriginalName"] = clusters_contigs.OriginalName

clusters["Large_enough"] = clusters.OriginalName.isin(big_bins)

del clusters_contigs

logging.info(f"Write reformated table to {output_culsters}")
clusters.to_csv(output_culsters, sep="\t")

clusters = clusters.query("Large_enough")

clusters["SampleBin"] = clusters.Sample + "_vamb_" + clusters.BinID

logging.info(f"Write cluster_attribution for samples")
for sample, cl in clusters.groupby("Sample"):
    sample_output_path = cluster_output_path.format(sample=sample)

    logging.debug(f"Write file {sample_output_path}")
    cl[["Contig", "SampleBin"]].to_csv(
        sample_output_path, sep="\t", index=False, header=False
    )


samples_without_bins = set(snakemake.params.samples).difference(set(clusters.Samples))

if len(samples_without_bins) > 0:
    logging.warning(
        "The following samples did't yield bins, I add longest contig to make the pipline continue:\n"
        + "\n".join(samples_without_bins)
    )

    for sample in samples_without_bins:
        sample_output_path = cluster_output_path.format(sample=sample)
        with open(sample_output_path, "w") as fout:
            fout.write(f"{sample}_0\t{sample}_vamb_1\n")
