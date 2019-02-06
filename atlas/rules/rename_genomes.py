
import os, shutil,sys
import argparse
from snakemake.io import glob_wildcards

from atlas import utils


def rename_genomes(input_folder,mapfile_genomes,mapfile_contigs,output_dir):

    file_name = f"{input_folder}/{{binid}}.fasta"
    bin_ids, = glob_wildcards(file_name)

    old2new_name= dict(zip(bin_ids,utils.gen_names_for_range(len(bin_ids),prefix='MAG')))
    os.makedirs(output_dir)


    with open(mapfile_contigs,'w') as out_contigs, open(mapfile_genomes,'w') as old2new_mapping_file :
        old2new_mapping_file.write(f"BinID\tMAG\n")
        for binid in bin_ids:

            fasta_in = file_name.format(binid=binid)
            new_name= old2new_name[binid]

            old2new_mapping_file.write(f"{binid}\t{new_name}\n")

            fasta_out = os.path.join(output_dir,f"{new_name}.fasta")

            # write names of contigs in mapping file
            with open(fasta_in) as ffi, open(fasta_out,'w') as ffo :
                Nseq=0
                for line in ffi:
                    if line[0]==">":
                        Nseq+=1
                        new_header=f'{new_name}_{Nseq}'
                        out_contigs.write(f"{new_header}\t{new_name}\n")
                        ffo.write(f">{new_header}\n")
                    else:
                        ffo.write(line)


if __name__ == "__main__":


    try:
        rename_genomes(
            input_folder=snakemake.input[0],
            output_dir=snakemake.output.dir,
            mapfile_genomes=snakemake.output.mapfile_genomes,
            mapfile_contigs=snakemake.output.mapfile_contigs
        )
    except NameError:
        p = argparse.ArgumentParser()
        p.add_argument("--input-folder")
        p.add_argument("--output-dir")
        p.add_argument("--mapfile-genomes")
        p.add_argument("--mapfile-contigs")
        args = vars(p.parse_args())
        rename_genomes(**args)
