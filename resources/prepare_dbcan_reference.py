import click
import re
import sys
from itertools import groupby


@click.command(context_settings={'help_option_names':['-h','--help']})
@click.argument("faminfo", click.Path(exists=True))
@click.argument("cazydb_fasta", click.Path(exists=True))
@click.argument("out_map")
@click.argument("out_fasta")
def prepare_ec(faminfo, cazydb_fasta, out_map, out_fasta):
    """

    faminfo:

        \b
        http://csbl.bmb.uga.edu/dbCAN/download/FamInfo.txt

    cazydb_fasta is:

        \b
        http://csbl.bmb.uga.edu/dbCAN/download/CAZyDB.07152016.fa

    """
    expazy = set()
    uniprot_to_uniparc = {}

    click.echo("parsing %s" % faminfo)
    fam_info = {}
    with open(faminfo) as fh:
        next(fh)
        for line in fh:
            toks = line.strip().split("\t")
            if "Unclassified" in toks[0]:
                toks[0] = toks[0].replace("Unclassified-", "")
            if toks[0] in fam_info:
                print(line)
                sys.exit("Non-unique family present in Family Info")
            # unsure how useful cazy_activities is here as it includes everything across all families
            fam_info[toks[0]] = toks[2]

    click.echo("parsing %s" % cazydb_fasta)
    cazy_fasta_map = {}
    with open(cazydb_fasta) as fh:
        for name, seq in read_fasta(fh):
            name_parts = name.split("|", 2)
            cazy_fasta_map[name_parts[0]] = {"seq":seq, "ecs":"" if len(name_parts) == 2 else name_parts[2], "cazy_family": name_parts[1]}

    with open(out_map, "w") as omap, open(out_fasta, "w") as ofa:
        for name, meta in cazy_fasta_map.items():
            # removes masked sequences
            if meta["seq"].islower(): continue
            print(format_fasta_record(name, meta["seq"]), file=ofa)
            # cazy_gene, cazy_family, cazy_class, cazy_ecs

            if not meta["cazy_family"] in fam_info:
                cazy_class = re.findall(r"([A-Z]+)", meta["cazy_family"])[0]
                if not cazy_class:
                    cazy_class = "NA"
            else:
                cazy_class = fam_info[meta["cazy_family"]]

            print(name, meta["cazy_family"], cazy_class, meta["ecs"], sep="\t", file=omap)


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
            name = line[1:].strip()
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
