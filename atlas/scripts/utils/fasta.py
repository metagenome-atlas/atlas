from Bio import SeqIO
from numpy import ceil
import os


def _make_test_fasta(test_file="test_ABC.fasta"):
    with open(test_file, "w") as f:
        for Number, Letter in enumerate("ATCG"):
            f.write(f">contig_{Number+1} description\n{Letter}\n")


def count_Nseq(fasta_file):
    """
    Counts number of sequences in a fasta file.
    >>> fasta_file='test_ABC.fasta'
    >>> _make_test_fasta(fasta_file) # makes fasta with a seq for each nucleotide
    >>> count_Nseq(fasta_file)
    4
    >>> os.remove(fasta_file)

    """
    i = 0
    with open(fasta_file) as f:
        for line in f:
            if line[0] == ">":
                i += 1
    return i


def split(fasta_file, maxSubsetSize, out_dir, simplify_headers=True):
    """
    Splits a fasta in subsets of size max maxSubsetSize.
    >>> fasta_file='test_ABC.fasta'
    >>> out_dir = 'test_outdit_doctest'
    >>> _make_test_fasta(fasta_file) # makes fasta with a seq for each nucleotide
    >>> split(fasta_file,3,out_dir,simplify_headers=True)
    >>> len(os.listdir(out_dir))
    2
    >>> count_Nseq('test_outdit_doctest/subset1.fasta')
    2
    >>> count_Nseq('test_outdit_doctest/subset2.fasta')
    2
    >>> split(fasta_file,3,out_dir,simplify_headers=True)
    Traceback (most recent call last):
        ...
    FileExistsError: [Errno 17] File exists: 'test_outdit_doctest'
    >>> import shutil; shutil.rmtree(out_dir)
    >>> os.remove(fasta_file)
    """

    N = count_Nseq(fasta_file)

    SubsetSize = int(ceil(N / ceil(N / maxSubsetSize)))
    extension = os.path.splitext(fasta_file)[-1]

    os.makedirs(out_dir)

    i, subset_n = 0, 0
    fout = None
    for i, seq in enumerate(SeqIO.parse(fasta_file, "fasta")):
        if (i % SubsetSize) == 0:
            subset_n += 1
            if fout is not None:
                fout.close()

            fout = open(f"{out_dir}/subset{subset_n}{extension}", "w")

        if simplify_headers:
            seq.description = ""
        SeqIO.write(seq, fout, "fasta")

    fout.close()


def header2origin(fasta_file, out, simplify_header=True):
    """
    Annotates a fasta file to it's filename:
    genome.fasta:
    >contig1 description
    ACTAC
    >contig2 description
    ACTAC
    ...

    becomes:
    contig1 genome
    contig2 genome

    input is a fasta filename:
    out is a filename or a stream

    """

    if type(out) == str:
        out_stream = open(out, "w")
    else:
        out_stream = out

    name = os.path.splitext(os.path.split(fasta_file)[-1])[0]

    # write names of contigs in mapping file
    with open(fasta_file) as f:
        for line in f:
            if line[0] == ">":
                header = line[1:].strip()
                if simplify_header:
                    header = header.split()[0]
                out_stream.write(f"{header}\t{name}\n")
    out_stream.flush()


if __name__ == "__main__":
    import doctest, shutil

    doctest.testmod()
