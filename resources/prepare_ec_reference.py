import click
import gzip
import sys
from collections import defaultdict
from itertools import groupby


@click.command(context_settings={'help_option_names':['-h','--help']})
@click.argument("enzyme_dat", click.Path(exists=True))
@click.argument("uniparc_map", click.Path(exists=True))
@click.argument("uniparc_fasta", click.Path(exists=True))
@click.argument("out_map")
@click.argument("out_fasta")
def prepare_ec(enzyme_dat, uniparc_map, uniparc_fasta, out_map, out_fasta):
    """

    enzyme_dat:

        \b
        ftp://ftp.expasy.org/databases/enzyme/enzyme.dat

    uniparc_map is the compressed version of was 'all' from:

        \b
        http://www.uniprot.org/uniparc/

    uniparc_fasta is:

        \b
        ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniparc/uniparc_active.fasta.gz

    """
    expazy = set()
    uniprot_to_uniparc = {}

    click.echo("parsing %s" % uniparc_map)
    with gzip.open(uniparc_map, 'rt') as fh:
        # Entry    Organisms    UniProtKB    First seen    Last seen    Length
        next(fh)
        for line in fh:
            toks = line.strip().split("\t")
            for uniprot_entry in toks[2].split(";"):
                uniprot_entry = uniprot_entry.strip()

                if not uniprot_entry or "obsolete" in uniprot_entry: continue
                if uniprot_entry in uniprot_to_uniparc:
                    print("wtf", uniprot_entry, line)
                    sys.exit(1)

                uniprot_to_uniparc[uniprot_entry] = toks[0]

    click.echo("parsing %s" % enzyme_dat)
    uniparc_mappings = defaultdict(list)
    with open(enzyme_dat) as fh:
        for key, group in groupby(fh, key=lambda x: x.startswith("//")):
            if not key:
                uniprot_entries = []
                recommended_name = "NA"
                ec_id = "NA"

                for line in group:
                    toks = line.strip().partition("   ")

                    if line.startswith("ID"):
                        ec_id = toks[2]

                    if line.startswith("DE"):
                        recommended_name = toks[2].strip(".")

                    if line.startswith("DR"):
                        for entry_and_name in toks[2].split(";"):
                            uniprot_entry = entry_and_name.split(",")[0].strip()

                            if uniprot_entry:
                                uniprot_entries.append(uniprot_entry)

                if ec_id and uniprot_entries:
                    for uniprot_entry in uniprot_entries:
                        uniparc_id = uniprot_to_uniparc[uniprot_entry]

                        uniparc_mappings[uniparc_id].append([uniprot_entry, ec_id, recommended_name])
                        expazy.add(uniparc_id)

    with open(out_map, "w") as ofh:
        for uniparc_id, meta in uniparc_mappings.items():
            uniprots = set()
            ecs = set()
            names = set()
            for name_list in meta:
                uniprots.add(name_list[0])
                ecs.add(name_list[1])
                names.add(name_list[2])
            print(uniparc_id, "|".join(uniprots), "|".join(ecs), "|".join(names), sep="\t", file=ofh)

    click.echo("parsing %s" % uniparc_fasta)
    with gzip.open(uniparc_fasta, 'rt') as fh, open(out_fasta, "w") as ofh:
        for name, seq in read_fasta(fh):
            if name in expazy:
                print(format_fasta_record(name, seq), file=ofh)


def read_fasta(fh):
    """Fasta iterator.

    Accepts file handle of .fasta and yields name and sequence.

    Args:
        fh (file): Open file handle of .fasta file

    Yields:
        tuple: name, sequence

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
    for header, group in groupby(fh, lambda line: line[0] == '>'):
        if header:
            line = next(group)
            name = line[1:].partition(" ")[0]
        else:
            seq = ''.join(line.strip() for line in group)
            yield name, seq


def format_fasta_record(name, seq, wrap=80):
    """Fasta __str__ method.

    Convert fasta name and sequence into wrapped fasta format.

    Args:
        name (str): name of the record
        seq (str): sequence of the record
        wrap (int): length of sequence per line

    Yields:
        tuple: name, sequence

    >>> format_fasta_record("seq1", "ACTG")
    ">seq1\nACTG"
    """
    record = ">" + name + "\n"
    if wrap:
        for i in range(0, len(seq), wrap):
            record += seq[i:i+wrap] + "\n"
    else:
        record += seq + "\n"
    return record.strip()


if __name__ == "__main__":
    prepare_ec()
