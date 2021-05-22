import contextlib
import gzip
import logging
import sqlite3
from itertools import groupby

from atlas import BLAST6
from atlas.blast import (
    BlastHits,
    Node,
    Tree,
    parse_blast_results_with_tree,
    process_orfs_with_tree,
)
from atlas.utils import gzopen


# TODO: update the docstring to reflect column changes
def refseq_parser(
    tsv,
    namemap,
    treefile,
    output,
    summary_method,
    aggregation_method,
    majority_threshold,
    min_identity,
    min_bitscore,
    min_length,
    max_evalue,
    max_hits,
    table_name,
    top_fraction,
):
    """Parse TSV (tabular BLAST output [-outfmt 6]), grabbing taxonomy metadata from ANNOTATION to
    compute LCAs.

    The BLAST hits are assumed to be sorted by query with decreasing bitscores (best alignment first):

        \b
        sort -k1,1 -k12,12rn tsv > sorted_tsv

    Annotation file should include your BLAST subject sequence ID, a function, a taxonomy name,
    the taxonomy ID, and the parent taxonomy ID. This file is generated from `prepare-refseq`:

        \b
        gi|507051347|ref|WP_016122340.1| two-component sensor histidine kinase Bacillus cereus 1396 86661
        gi|507052147|ref|WP_016123132.1| two-component sensor histidine kinase Bacillus cereus 1396 86661
        gi|507053266|ref|WP_016124222.1| two-component sensor histidine kinase Bacillus cereus 1396 86661

    The RefSeq function is always derived from the best BLAST hit.

    The output will give contig, ORF ID, the lineage assigned to the contig based on
    --aggregation-method, the probability of error (erfc), taxonomy assigned to the ORF, the
    best hit's product, the best hit's evalue, and the best hit's bitscore:

        \b
        contig     orf          taxonomy erfc orf_taxonomy refseq               refseq_evalue refseq_bitscore
        k121_52126 k121_52126_1 root     1.0  root         hypothetical protein 1.0e-41       176.8

    Args:
        tsv (str): blast hits file path in default tabular format
        namemap (str): sqlite database file path
        treefile (str): file path to NCBI tree file
        output (str): :py:class:`click.File` in write mode
        summary_method (str): either 'majority' or 'best'; summary method for annotating ORFs; when majority and there is no majority, best is used
        aggregation_method (str): 'lca', 'lca-majority', or 'majority; summary method for aggregating ORF taxonomic assignments to contig level assignment; 'lca' will result in most stringent, least specific assignments
        majority_threshold (float): constitutes a majority fraction at tree node for 'lca-majority' ORF aggregation method
        min_identity (int): minimum allowable percent ID of BLAST hit
        min_bitscore (int): minimum allowable bitscore of BLAST hit; 0 disables
        min_length (int): minimum allowable BLAST alignment length
        max_evalue (float): maximum allowable e-value of BLAST hit
        top_fraction (float): filters ORF BLAST hits before finding majority by only keep hits within this fraction, e.g. 0.98, of the highest bitscore
        max_hits (int): maximum number of BLAST hits to consider when summarizing ORFs as a majority
        table_name (str): table name within namemap database; expected columns are listed above

    """
    logging.info("Parsing %s" % tsv)
    tree = Tree(treefile)
    orf_assignments = parse_blast_results_with_tree(
        tsv,
        namemap,
        summary_method=summary_method,
        tree=tree,
        min_identity=min_identity,
        min_bitscore=min_bitscore,
        min_length=min_length,
        max_evalue=max_evalue,
        max_hits_per_orf=max_hits,
        table_name=table_name,
        top_fraction_of_hits=top_fraction,
    )
    logging.info("Assigning taxonomies to contigs using %s" % aggregation_method)
    process_orfs_with_tree(
        orf_assignments,
        tree,
        output,
        aggregation_method,
        majority_threshold,
        table_name,
    )
    logging.info("Complete")


def read_fasta(fh):
    """Fasta iterator.

    Accepts file handle of .fasta and yields name and sequence.

    Args:
        fh (file): Open file handle of .fasta file

    Yields:
        tuple: name, sequence

    Example:
        >>> import os
        >>> from itertools import groupby
        >>> f = open("test.fasta", 'w')
        >>> f.write("@seq1\nACTG")
        >>> f.close()
        >>> f = open("test.fastq")
        >>> for name, seq in read_fastq(f):
                assert name == "seq1"
                assert seq == "ACTG"
        >>> f.close()
        >>> os.remove("test.fasta")

    """
    for header, group in groupby(fh, lambda line: line[0] == ">"):
        if header:
            line = next(group)
            name = line[1:].strip()
        else:
            seq = "".join(line.strip() for line in group)
            yield name, seq
