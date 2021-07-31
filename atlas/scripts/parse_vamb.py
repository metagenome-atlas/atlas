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

vamb_cluster_file= "Crossbinning/vamb/clustering/clusters.tsv"# HACK:
output_culsters = "Crossbinning/vamb/Clusters.tsv"


logging.info(f"Load vamb cluster file {vamb_cluster_file}")
clusters = pd.read_table(vamb_cluster_file, header=None)

clusters.columns=['OriginalName', 'Contig']


clusters['Sample'] = clusters.Contig.str.rsplit('_', n=1, expand=True)[0]

# Calculate cluster size drop singletons
Cluster_Size=  clusters.OriginalName.value_counts()
Cluster_Size= Cluster_Size.loc[lambda x: >1]
logging.info(Cluster_Size)

# Create Cluster names
unique_clusters =Cluster_Size.index
rename_bins= dict(zip( unique_clusters , gen_names_for_range(N= len(unique_clusters), prefix= "Cluster" )))
clusters['Cluster'] = clusters.OriginalName.map( rename_bins)

# drop clusters with no new name = size below
clusters.dropna(inplace=True)

# create unique bins for each sample
clusters['SampleBin'] = "Vamb_" +  clusters.Cluster + "_" + clusters.Sample

clusters.set_index('Contigs',drop=True, inplace=True)

logging.info(f"Write reformated table to {output_culsters}")
clusters.to_csv(output_culsters,sep='\t')

logging.info(f"Write cluster_attribution for samples")
#output_tables = dict(zip(snakemake.params.Samples, snakemake.output.cluster_attributions))
output_tables={} # HACK:
for sample, cl in clusters.groupby('Sample'):

    output_tables[sample] = f"{sample}/binning/Vamb/cluster_attribution.tsv"# HACK:

    logging.debug(f"Write file {output_tables[sample]}")
    cl.SampleBin.to_csv(output_tables[sample], sep='\t', header=False )
