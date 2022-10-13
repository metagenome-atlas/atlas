from multiprocessing import Pool
import pandas as pd
import os, sys
from .io import simplify_path, simply_open
from itertools import groupby
import numpy as np
import gzip as gz

def verify_dna(sequence, is_upper):

    if not is_upper:
        sequence=sequence.upper()
    letters_used =set(sequence)

    alphabet= set(list("ATCGN"))

    additional_letters= letters_used -alphabet

    if len(additional_letters)>0:
        raise Exception(f"Sequence contains additional letters that are not DNA {additional_letters}")

def get_stats_from_lengths(lengths):

    sorted_lengths = sorted(lengths, reverse=True)
    csum = np.cumsum(sorted_lengths)

    Total_length = int(sum(lengths))
    N = len(lengths)

    n2 = int(Total_length / 2)

    # get index for cumsum >= N/2
    csumn2 = min(csum[csum >= n2])
    ind = int(np.where(csum == csumn2)[0][0])

    N50 = sorted_lengths[ind]

    return Total_length, N, N50


def genome_stats(fasta_file, number_of_n_for_split=10):
    """Get genome stats from a fasta file. Outputs a tuple with:
       name,Length, n_seq,N50
    """

    try:

        name = simplify_path(fasta_file)

        scaffold_lengths = []
        contig_lengths = []
        ambigious_bases = 0

        with simply_open(fasta_file, "r") as fasta:
            ## parse each sequence by header: groupby(data, key)
            faiter = (x[1] for x in groupby(fasta, lambda line: line[0] == ">"))

            for record in faiter:
                # reccord contains header
                ## join sequence lines
                sequence = ''.join(s.strip() for s in faiter.__next__())
                sequence = sequence.upper()

                verify_dna(sequence, is_upper=True)

                # count ambigous bases
                ambigious_bases+= sequence.count("N")


                # get set of scaffold lengths
                scaffold_lengths.append(len(sequence))

                # split on Ns and get set of contig lengths.
                contig_lengths += [
                    len(contig) for contig in sequence.split("N"*number_of_n_for_split)
                ]

        Length_scaffolds, N_scaffolds, N50 = get_stats_from_lengths(scaffold_lengths)

        Length_contigs, N_contigs, _ = get_stats_from_lengths(contig_lengths)

    except Exception as e:
        raise Exception(f"Error in calculating stats of {fasta_file}") from e

    return {"File":name, 
            "Length_scaffolds":Length_scaffolds,
            "N_scaffolds":N_scaffolds, 
            "N50":N50, 
            "Length_contigs":Length_contigs,
            "N_contigs": N_contigs,
            "Ambigious_bases": ambigious_bases}



def get_many_genome_stats(filenames, output_filename, threads=1):
    """Small function to calculate total genome length and N50
    """

    pool = Pool(threads)

    results = pool.map(genome_stats, filenames)
    Stats = pd.DataFrame(results).rename({"Length_scaffolds":"Length"})

    Stats.to_csv(output_filename, sep="\t", index=False)
