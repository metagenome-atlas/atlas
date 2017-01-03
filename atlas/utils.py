import logging
import os
from collections import Counter
from math import log, erfc, sqrt
from snakemake.io import load_configfile


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
                continue
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
