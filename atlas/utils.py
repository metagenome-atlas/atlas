import logging
import os
import re
import sys
from collections import Counter, OrderedDict
from itertools import groupby
from math import log, erfc, sqrt
from snakemake.io import load_configfile


gzopen = lambda f: gzip.open(f, mode="rt") if f.endswith(".gz") else open(f)


def validate_assembly_config(config):
    c = load_configfile(config)
    valid = True
    if not "samples" in c:
        logging.critical("'samples' is not defined in %s" % config)
        valid = False
    try:
        if len(c["samples"].keys()) == 0:
            logging.critical("no samples are defined under 'samples' in %s" % config)
            valid = False
        for sample, meta in c["samples"].items():
            if sample == "coassemblies":
                for coassembly, file_list in sample["coassemblies"].items():
                    for sample in file_list:
                        if sample not in c["samples"]:
                            logging.critical("Sample %s under coassembly %s is not a defined sample in the configuration" % (sample, coassembly))
                            valid = False
            if not "path" in meta:
                logging.critical("'path' is not set for sample %s" % sample)
                valid = False
                continue
            for f in meta["path"]:
                if not os.path.exists(f):
                    logging.critical("%s does not exist for sample %s" % (f, sample))
                    valid = False
    except KeyError:
        pass

    if not "preprocessing" in c:
        logging.critical("'preprocessing' is not defined in %s" % config)
        valid = False
    try:
        if not "adapters" in c["preprocessing"]:
            logging.critical("'adapters' is not defined under 'preprocessing' in %s" % config)
            valid = False
        for f in c["preprocessing"]["adapters"].split(","):
            if not os.path.exists(f):
                logging.critical("adapters file [%s] does not exist" % f)
                valid = False
        if not "contamination" in c["preprocessing"]:
            logging.critical("'contamination' is not defined under 'preprocessing' in %s" % config)
            valid = False
        if not "references" in c["preprocessing"]["contamination"]:
            logging.critical("'references' is not defined under 'contamination' in %s" % config)
            valid = False
        if not "rRNA" in c["preprocessing"]["contamination"]["references"]:
            logging.critical("'rRNA' is not a defined contamination reference in %s" % config)
        for ref, f in c["preprocessing"]["contamination"]["references"].items():
            if not os.path.exists(f):
                logging.critical("contamination reference file [%s] does not exist for %s" % (f, ref))
                valid = False
        if not "normalization" in c["preprocessing"]:
            logging.critical("'normalization' is not defined in %s" % config)
            valid = False
    except KeyError:
        pass

    if not "assembly" in c:
        logging.critical("'assembly' is not defined in %s" % config)
        valid = False
    try:
        if not c["assembly"]["assembler"] == "megahit" and not c["assembly"]["assembler"] == "spades":
            logging.critical("'assembler' entry [%s] is not a supported assembler" % c["assembly"]["assembler"])
            valid = False
    except KeyError:
        pass

    if not "annotation" in c:
        logging.critical("'annotation' is not defined in %s" % config)
        valid = False
    try:
        if not "references" in c["annotation"]:
            logging.critical("'references' is not defined in %s" % config)
            valid = False
        # TODO: decide which if any of these are required
        if not "eggnog" in c["annotation"]["references"]:
            logging.critical("'eggnog' is not defined as an annotation reference")
            valid = False
        if not "refseq" in c["annotation"]["references"]:
            logging.critical("'refseq' is not defined as an annotation reference")
            valid = False
        if not "cazy" in c["annotation"]["references"]:
            logging.critical("'cazy' is not defined as an annotation reference")
            valid = False
        if not "expazy" in c["annotation"]["references"]:
            logging.critical("'expazy' is not defined as an annotation reference")
            valid = False
        if not os.path.exists(c["annotation"]["references"]["eggnog"]["namemap"]):
            logging.critical("namemap reference file [%s] does not exist" % c["annotation"]["references"]["eggnog"]["namemap"])
            valid = False
        if not os.path.exists(c["annotation"]["references"]["eggnog"]["dmnd"]):
            logging.critical("fasta reference file [%s] does not exist" % c["annotation"]["references"]["eggnog"]["fasta"])
            valid = False
        if not os.path.exists(c["annotation"]["references"]["refseq"]["namemap"]):
            logging.critical("namemap reference file [%s] does not exist" % c["annotation"]["references"]["refseq"]["namemap"])
            valid = False
        if not os.path.exists(c["annotation"]["references"]["refseq"]["tree"]):
            logging.critical("tree reference file [%s] does not exist" % c["annotation"]["references"]["refseq"]["tree"])
            valid = False
        if not os.path.exists(c["annotation"]["references"]["refseq"]["dmnd"]):
            logging.critical("fasta reference file [%s] does not exist" % c["annotation"]["references"]["refseq"]["fasta"])
            valid = False
        if not os.path.exists(c["annotation"]["references"]["cazy"]["namemap"]):
            logging.critical("namemap reference file [%s] does not exist" % c["annotation"]["references"]["cazy"]["namemap"])
            valid = False
        if not os.path.exists(c["annotation"]["references"]["cazy"]["dmnd"]):
            logging.critical("fasta reference file [%s] does not exist" % c["annotation"]["references"]["cazy"]["fasta"])
            valid = False
        if not os.path.exists(c["annotation"]["references"]["expazy"]["namemap"]):
            logging.critical("namemap reference file [%s] does not exist" % c["annotation"]["references"]["expazy"]["namemap"])
            valid = False
        if not os.path.exists(c["annotation"]["references"]["expazy"]["dmnd"]):
            logging.critical("fasta reference file [%s] does not exist" % c["annotation"]["references"]["expazy"]["fasta"])
            valid = False
    except KeyError:
        pass

    return valid


def index_of_list_items(lists):
    """Find the highest index for items among list of lists and return ordered dictionary sorted by
    decreasing index value.

    Args:
        lists (list): list of lists

    Returns:
        OrderedDict of list items, sorted by depth with deepest list item first

    Example:
        >>> lineages = [['1', '131567', '99999999', '2'],
                        ['1', '131567', '99999999', '2', '1224', '1236', '72274', '135621', '286'],
                        ['1', '131567', '99999999', '2', '1224'],
                        ['1', '131567', '99999999', '2', '1224', '1236']]
        >>> index_of_list_items(lineages)
        OrderedDict([('286', 8),
                     ('135621', 7),
                     ('72274', 6),
                     ('1236', 5),
                     ('1224', 4),
                     ('2', 3),
                     ('99999999', 2),
                     ('131567', 1),
                     ('1', 0)])

    """
    indexes = {}
    for l in lists:
        for i, item in enumerate(l):
            if item in indexes:
                if indexes[item] < i:
                    indexes[item] = i
            else:
                indexes[item] = i
    return OrderedDict(sorted(indexes.items(), key=lambda t: t[1], reverse=True))


def nettleton_pvalue(items, key):
    """Calculate pvalue based on Nettleton result.

    Args:
        items (list): list of items
        key (string): key within items you're testing

    Returns:
        float

    Raises:
        AssertionError when key is not present in items
    """
    assert key in items
    if len(items) == 1:
        return 1

    item_counts = Counter(items)
    if item_counts.most_common()[0][0] == key:
        try:
            # take second entries count value
            max_count = item_counts.most_common()[1][1]
        except IndexError:
            max_count = 0
    else:
        # take first entries count value
        max_count = item_counts.most_common()[0][1]
    if item_counts[key] <= max_count:
        return 1
    else:
        try:
            t = 2 * (max_count * log((2 * max_count) / (max_count + item_counts[key])) \
                  + (item_counts[key] * log((2 * item_counts[key] / \
                                                         (max_count + item_counts[key])))))
        except ValueError:
            # max_count = 0
            t = 2 * (item_counts[key] * log((2 * item_counts[key] / \
                                                         (max_count + item_counts[key]))))
        return erfc(sqrt(t / 2))


def gff_to_gtf(gff_in, gtf_out):
    t = re.compile(r'ID=[0-9]+_([0-9]+);')
    with open(gtf_out, "w") as fh, open(gff_in) as gff:
        for line in gff:
            if line.startswith("#"): continue
            toks = line.strip().split("\t")
            orf = t.findall(toks[-1])[0]
            gene_id = toks[0] + "_" + orf
            toks[-1] = 'gene_id "%s"; %s' % (gene_id, toks[-1])
            print(*toks, sep="\t", file=fh)


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


def print_fasta_record(name, seq, out_handle=sys.stdout, wrap=80):
    """Print a fasta record accounting for line wraps in the sequence.

    Args:
        name (str): name or header for fasta entry
        seq (str): sequence
        out_handle (Optional): open file handle in which to write or stdout
        wrap (Optional[int]) : line width of fasta sequence; None is supported for no wrapping

    """
    print('>', name, sep='', file=out_handle)
    if wrap:
        for i in range(0, len(seq), wrap):
            print(seq[i:i + wrap], file=out_handle)
    else:
        print(seq, file=out_handle)


def split_fasta(fasta, chunk_size=250000):
    chunk_size = int(chunk_size)
    fasta = os.path.expanduser(fasta)
    root, ext = os.path.splitext(fasta)

    file_idx = 0
    with open(fasta) as f:
        for i, (name, seq) in enumerate(read_fasta(f)):
            if i % chunk_size == 0:
                if i == 0:
                    ofh = open("{root}_{idx}{ext}".format(root=root, idx=file_idx, ext=ext), "w")
                    print_fasta_record(name, seq, out_handle=ofh)
                else:
                    ofh.close()
                    file_idx += 1
                    ofh = open("{root}_{idx}{ext}".format(root=root, idx=file_idx, ext=ext), "w")
                    print_fasta_record(name, seq, out_handle=ofh)
            else:
                print_fasta_record(name, seq, out_handle=ofh)
    ofh.close()
