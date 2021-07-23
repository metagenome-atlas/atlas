import os, shutil, sys
import argparse
from snakemake.io import glob_wildcards

from atlas import utils
import pandas as pd


def rename_genomes(
    input_folder, mapfile_genomes, mapfile_contigs, output_dir, rename_contigs=True
):

    file_name = f"{input_folder}/{{binid}}.fasta"
    (bin_ids,) = glob_wildcards(file_name)

    old2new_name = dict(
        zip(bin_ids, utils.gen_names_for_range(len(bin_ids), prefix="MAG"))
    )
    os.makedirs(output_dir)

    with open(mapfile_contigs, "w") as out_contigs, open(
        mapfile_genomes, "w"
    ) as old2new_mapping_file:
        old2new_mapping_file.write(f"BinID\tMAG\n")
        for binid in bin_ids:

            fasta_in = file_name.format(binid=binid)
            new_name = old2new_name[binid]

            old2new_mapping_file.write(f"{binid}\t{new_name}\n")

            fasta_out = os.path.join(output_dir, f"{new_name}.fasta")

            # write names of contigs in mapping file
            with open(fasta_in) as ffi, open(fasta_out, "w") as ffo:
                Nseq = 0
                for line in ffi:
                    if line[0] == ">":
                        Nseq += 1
                        if rename_contigs:
                            new_header = f"{new_name}_{Nseq}"
                        else:
                            new_header = line[1:].strip().split()[0]
                        out_contigs.write(f"{new_header}\t{new_name}\n")
                        ffo.write(f">{new_header}\n")
                    else:
                        ffo.write(line)


def genome2cluster(Drep_folder):

    Cdb = pd.read_csv(os.path.join(Drep_folder, "..", "data_tables", "Cdb.csv"))
    Cdb.index = Cdb.genome  # location changes

    Wdb = pd.read_csv(os.path.join(Drep_folder, "..", "data_tables", "Wdb.csv"))
    Wdb.index = Wdb.cluster
    map_genome2cluster = Cdb.secondary_cluster.map(Wdb.genome)

    return map_genome2cluster


def get_mapfile_bins(mapfile_bins, dereplication, old2new):

    Genome_map = genome2cluster(dereplication)

    Genome_map.index = Genome_map.index.str.replace(".fasta", "")
    Genome_map = Genome_map.str.replace(".fasta", "")

    old2new_name = pd.read_csv(old2new, index_col=0, squeeze=True, sep="\t")
    Genome_map = Genome_map.map(old2new_name)

    Genome_map.sort_values(inplace=True)
    Genome_map.name = "MAG"

    Genome_map.to_csv(mapfile_bins, sep="\t", header=True)


if __name__ == "__main__":

    try:
        log = open(snakemake.log[0], "w")
        sys.stderr = log
        sys.stdout = log

        rename_genomes(
            input_folder=snakemake.input.genomes,
            output_dir=snakemake.output.dir,
            mapfile_genomes=snakemake.output.mapfile_genomes,
            mapfile_contigs=snakemake.output.mapfile_contigs,
            rename_contigs=snakemake.params.rename_contigs,
        )

        get_mapfile_bins(
            mapfile_bins=snakemake.output.mapfile_bins,
            dereplication=snakemake.input.genomes,
            old2new=snakemake.output.mapfile_genomes,
        )

    except NameError:
        p = argparse.ArgumentParser()
        p.add_argument("--input-folder")
        p.add_argument("--output-dir")
        p.add_argument("--mapfile-genomes")
        p.add_argument("--mapfile-contigs")
        args = vars(p.parse_args())
        rename_genomes(**args)
