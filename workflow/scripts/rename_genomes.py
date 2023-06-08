#! /usr/bin/env python

import sys, os
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


# start


from snakemake.io import glob_wildcards

from atlas import utils
import pandas as pd

# Bin     Filename        Proteins
paths = pd.read_csv(snakemake.input.paths, sep="\t", index_col=0).Filename
# genome  SpeciesNr       Species Representative
mapping = pd.read_csv(
    snakemake.input.mapping_file,
    sep="\t",
    index_col=0,
).squeeze()


# standardize names of representatives
# MAG001 ....
representatives = mapping.Representative.unique()
old2new_name = dict(
    zip(representatives, utils.gen_names_for_range(len(representatives), prefix="MAG"))
)
mapping["MAG"] = mapping.Representative.map(old2new_name)


# write cluster attribution
mapping[["MAG", "Representative"]].to_csv(
    snakemake.output.mapfile_allbins2mag, sep="\t", header=True
)

# write out old2new ids
old2new = mapping.loc[representatives, "MAG"]
old2new.index.name = "Representative"
old2new.to_csv(snakemake.output.mapfile_old2mag, sep="\t", header=True)


#### Write genomes and contig to genome mapping file
output_dir = snakemake.output.dir
mapfile_contigs = snakemake.output.mapfile_contigs
rename_contigs = snakemake.params.rename_contigs


os.makedirs(output_dir)

with open(mapfile_contigs, "w") as out_contigs:
    for rep in representatives:
        fasta_in = paths.loc[rep]
        new_name = old2new.loc[rep]

        fasta_out = os.path.join(output_dir, f"{new_name}.fasta")

        # write names of contigs in mapping file
        with open(fasta_in) as ffi, open(fasta_out, "w") as ffo:
            Nseq = 0
            for line in ffi:
                # if header line
                if line[0] == ">":
                    Nseq += 1

                    if rename_contigs:
                        new_header = f"{new_name}_{Nseq}"
                    else:
                        new_header = line[1:].strip().split()[0]

                    # write to contig to mapping file
                    out_contigs.write(f"{new_header}\t{new_name}\n")
                    # write to fasta file
                    ffo.write(f">{new_header}\n")
                else:
                    ffo.write(line)


# rename quality
def rename_quality(quality_in, quality_out, old2new_name):
    Q = pd.read_csv(quality_in, index_col=0, sep="\t")

    Q = Q.loc[old2new_name.keys()].rename(index=old2new_name)

    Q.to_csv(quality_out, sep="\t")


rename_quality(
    quality_in=snakemake.input.genome_info,
    quality_out=snakemake.output.genome_info,
    old2new_name=old2new_name,
)
