#!/usr/bin/env python
# coding=utf-8
"""
Assumes BLAST output (-outfmt 6) is sorted by query then best hit (the default behavior). If
you're not sure, you can run something like:

sort -k1,1 -k12,12nr blast_output.tsv > blast_output.sorted.tsv

Metat_1000000_contigs_6828_1

Metat_1000000_contigs_6828_1	gi|495355777|ref|WP_008080498.1|	96	340	12	0	0	340	442	782	3.22e-209	1718
Metat_1000000_contigs_6828_1	gi|499651245|ref|WP_011331979.1|	99	340	2	0	0	340	442	782	5.71e-215	1763
Metat_1000000_contigs_6828_1	gi|495227881|ref|WP_007952652.1|	100	340	0	0	0	340	442	782	7.28e-216	1770

brow015@we30325:~/devel/blastlca
$ grep "gi|495227881|ref|WP_007952652.1|" refseq-nr-2014-01-18-names.txt
>gi|495227881|ref|WP_007952652.1| molybdopterin-dependent oxidoreductase alpha subunit [Pseudomonas sp. GM25]

brow015@we30325:~/devel/blastlca
$ grep "gi|499651245|ref|WP_011331979.1|" refseq-nr-2014-01-18-names.txt
>gi|499651245|ref|WP_011331979.1| oxidoreductase alpha (molybdopterin) subunit [Pseudomonas fluorescens]





oldquery = [["Pseudomonas sp. G5(2012)"],["Pseudomonas sp. HPB0071"],["Azotobacter vinelandii"],["Methylohalobius crimeensis"],["Pelobacter propionicus"],["Candidatus Methylomirabilis"],["Nevskia ramosa"],["Methylopila sp. M107"],["delta proteobacterium NaphS2"],["Hydrocarboniphaga effusa"]]
oldtree.getTaxonomy(oldquery)
'not Bacteria Haeckel 1894'

newtree = Tree("/Users/brow015/devel/blastlca/ncbi_taxonomy_tree.txt")
dd = defaultdict(lambda: defaultdict(BlastHit))
dd["Metat_1000000_contigs_60988_0"]["0"].add("Pseudomonas sp. G5(2012)", "1225")
dd["Metat_1000000_contigs_60988_0"]["0"].add("Pseudomonas sp. HPB0071", "889")
dd["Metat_1000000_contigs_60988_0"]["0"].add("Azotobacter vinelandii", 887)
dd["Metat_1000000_contigs_60988_0"]["0"].add("Methylohalobius crimeensis", 671)
dd["Metat_1000000_contigs_60988_0"]["0"].add("Pelobacter propionicus", 668)
dd["Metat_1000000_contigs_60988_0"]["0"].add("Candidatus Methylomirabilis", 638)
dd["Metat_1000000_contigs_60988_0"]["0"].add("Nevskia ramosa", 631)
dd["Metat_1000000_contigs_60988_0"]["0"].add("Methylopila sp. M107", 600)
dd["Metat_1000000_contigs_60988_0"]["0"].add("delta proteobacterium NaphS2", 598)
dd["Metat_1000000_contigs_60988_0"]["0"].add("Hydrocarboniphaga effusa", 593)
dd["Metat_1000000_contigs_60988_0"]["0"].add("TEST", 500)
dd['Metat_1000000_contigs_60988_0']['0'].taxonomies.reverse()
newtree.lca(dd['Metat_1000000_contigs_60988_0']['0'].taxonomies)

dd['Metat_1000000_contigs_60988_0']['0'].lca(newtree)


newtree = Tree("/Users/brow015/devel/blastlca/ncbi_taxonomy_tree.txt")
newtree.lca_star(['purple photosynthetic bacteria and relatives', 'gamma subgroup', 'not Bacteria Haeckel 1894', 'RNA similarity group I'])



contigs['Metat_1000000_contigs_60988']['4'].add('Pseudomonas sp. HPB0071', 1431, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Pseudomonas sp. HPB0071', 1377, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Pseudomonas mandelii', 1333, 10, 1)
contigs['Metat_1000000_contigs_60988']['6'].add('Pseudomonas mandelii', 1180, 10, 1)
contigs['Metat_1000000_contigs_60988']['6'].add('Pseudomonas sp. HPB0071', 1167, 10, 1)
contigs['Metat_1000000_contigs_60988']['0'].add('Pseudomonas sp. HPB0071', 889, 10, 1)
contigs['Metat_1000000_contigs_60988']['1'].add('Pseudomonas sp. Chol1', 844, 10, 1)
contigs['Metat_1000000_contigs_60988']['1'].add('Pseudomonas sp. HPB0071', 809, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Pseudomonas agarici', 1326, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Pseudomonas veronii', 1290, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Serratia marcescens', 1283, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Pantoea sp. A4', 1253, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia coli', 1237, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Shigella dysenteriae', 1234, 10, 1)
contigs['Metat_1000000_contigs_60988']['6'].add('Pseudomonas agarici', 1175, 10, 1)
contigs['Metat_1000000_contigs_60988']['6'].add('Pseudomonas veronii', 1111, 10, 1)
contigs['Metat_1000000_contigs_60988']['1'].add('Pseudomonas sp. 313', 799, 10, 1)
contigs['Metat_1000000_contigs_60988']['4'].add('Pseudomonas luteola', 1420, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Serratia marcescens', 1277, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia coli', 1238, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia coli', 1238, 10, 1)
contigs['Metat_1000000_contigs_60988']['6'].add('Pseudomonas luteola', 1114, 10, 1)
contigs['Metat_1000000_contigs_60988']['1'].add('Pseudomonas stutzeri', 856, 10, 1)
contigs['Metat_1000000_contigs_60988']['1'].add('Pseudomonas luteola', 809, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Salmonella enterica', 272, 10, 1)
contigs['Metat_1000000_contigs_60988']['4'].add('Pseudomonas sp. G5(2012)', 1939, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Pseudomonas sp. G5(2012)', 1900, 10, 1)
contigs['Metat_1000000_contigs_60988']['5'].add('Pseudomonas sp. G5(2012)', 1716, 10, 1)
contigs['Metat_1000000_contigs_60988']['6'].add('Pseudomonas sp. G5(2012)', 1661, 10, 1)
contigs['Metat_1000000_contigs_60988']['2'].add('Pseudomonas sp. G5(2012)', 1421, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Serratia liquefaciens', 1279, 10, 1)
contigs['Metat_1000000_contigs_60988']['1'].add('Pseudomonas sp. G5(2012)', 1228, 10, 1)
contigs['Metat_1000000_contigs_60988']['0'].add('Pseudomonas sp. G5(2012)', 1225, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia albertii', 1217, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia albertii', 1217, 10, 1)
contigs['Metat_1000000_contigs_60988']['7'].add('Pseudomonas sp. G5(2012)', 1150, 10, 1)
contigs['Metat_1000000_contigs_60988']['1'].add('Pseudomonas stutzeri', 878, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Serratia marcescens', 1283, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Serratia marcescens', 1273, 10, 1)
contigs['Metat_1000000_contigs_60988']['1'].add('Pseudomonas stutzeri', 816, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia coli', 1237, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia sp. TW11588', 1234, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia coli', 1234, 10, 1)
contigs['Metat_1000000_contigs_60988']['1'].add('Pseudomonas stutzeri subgroup', 894, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia coli', 618, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Enterobacteriaceae', 246, 10, 1)
contigs['Metat_1000000_contigs_60988']['4'].add('Azotobacter vinelandii', 1462, 10, 1)
contigs['Metat_1000000_contigs_60988']['2'].add('Azotobacter vinelandii', 992, 10, 1)
contigs['Metat_1000000_contigs_60988']['1'].add('Pseudomonas stutzeri', 892, 10, 1)
contigs['Metat_1000000_contigs_60988']['1'].add('Pseudomonas stutzeri', 879, 10, 1)
contigs['Metat_1000000_contigs_60988']['1'].add('Pseudomonas stutzeri', 878, 10, 1)
contigs['Metat_1000000_contigs_60988']['7'].add('Azotobacter vinelandii', 780, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Azotobacter vinelandii', 1355, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Serratia marcescens', 1284, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia coli', 1240, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia coli', 1237, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia coli', 1237, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia coli', 1233, 10, 1)
contigs['Metat_1000000_contigs_60988']['6'].add('Azotobacter vinelandii', 1132, 10, 1)
contigs['Metat_1000000_contigs_60988']['1'].add('Pseudomonas stutzeri', 894, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Shigella flexneri', 403, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Erwinia billingiae', 1248, 10, 1)
contigs['Metat_1000000_contigs_60988']['0'].add('Azotobacter vinelandii', 887, 10, 1)
contigs['Metat_1000000_contigs_60988']['1'].add('Azotobacter vinelandii', 849, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Enterobacteriaceae bacterium 9_2_54FAA', 1252, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Hafnia alvei', 1250, 10, 1)
contigs['Metat_1000000_contigs_60988']['1'].add('Pseudomonas stutzeri', 881, 10, 1)
contigs['Metat_1000000_contigs_60988']['1'].add('Pseudomonas stutzeri', 877, 10, 1)
contigs['Metat_1000000_contigs_60988']['1'].add('Pseudomonas stutzeri', 850, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia coli', 658, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Shimwellia blattae', 1240, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Serratia sp. M24T3', 1239, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia coli', 1234, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia coli', 1233, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia coli', 1233, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia coli', 1232, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia coli', 1231, 10, 1)
contigs['Metat_1000000_contigs_60988']['8'].add('Escherichia coli', 1140, 10, 1)
contigs['Metat_1000000_contigs_60988']['5'].add('Shigella boydii', 373, 10, 1)


"""
import argparse
import bisect
import logging
import math
import os
import re
from collections import Counter, defaultdict, deque, OrderedDict
from math import log, sqrt, erfc


class Tree(object):

    def __init__(self, tax_tree):
        """
        Raises:
            AssertionError when child and parent IDs are equal when not root
        """
        self.tree = defaultdict(dict)

        with open(tax_tree) as tree_fh:
            for line in tree_fh:
                toks = line.strip().split("\t")
                if not toks[1] == '1' and not toks[2] == '1':
                    assert not toks[1] == toks[2]
                if not len(toks) == 3:
                    logging.warning("Line [%s] does not have NAME, ID, PARENTID" % line.strip())
                    continue
                self.add_node(toks[0], toks[1], toks[2])

    def add_node(self, taxonomy, node_id, parent_id):
        # taxonomy string to node mapping
        self.tree[taxonomy] = Node(taxonomy, node_id, parent_id)
        # taxonomy id to node mapping
        self.tree[node_id] = self.tree[taxonomy]

    def taxonomy_id_to_name(self, id_map, key):
        """
        Args:
            id_map (dict): the keys are taxonomy IDs with values of names you'd like to use instead
                of the names in the taxonomy tree file
            key (str): taxonomy ID used to search `id_map`

        Returns:
            string of mapped key along with key, e.g. New Name (1224), or if they key is missing
                from `id_map`, it is translated into taxonomy name and used instead
        """
        val = "{name} ({digit})"
        if key in id_map:
            return val.format(name=id_map[key], digit=key)
        else:
            return val.format(name=self.tree[key].taxonomy, digit=key)

    def lca(self, taxonomies, threshold=1.):
        """Returns the taxonomy of the LCA and optionally only use the top fraction of hits.

        Args:
            taxonomies: list of taxonomies (ID or name); when using threshold < 1 they should be
                ordered by decreasing bitscore
            threshold (Optional[float]): 0-1; threshold fraction of hits to be factored into lca

        Returns:
            str: taxonomy of LCA
        """
        if threshold > 1:
            threshold = 1
        elif threshold < 0.01:
            # 10% as the minimum
            threshold = 0.1

        count_target = len(taxonomies) * threshold
        count_taxonomies = defaultdict(int)

        for taxonomy in taxonomies:
            try:
                current_taxonomy = self.tree[taxonomy].node_id
            except AttributeError:
                # dict when key not present
                # taxonomy represented in the reference database, but is not present in the tree
                continue
            while not current_taxonomy == "1":
                count_taxonomies[current_taxonomy] += 1
                if count_taxonomies[current_taxonomy] >= count_target:
                    return self.tree[current_taxonomy].taxonomy
                # traverse up tree
                current_taxonomy = self.tree[current_taxonomy].parent_id
        return "root"

    def filter_taxonomy_list(self, taxonomy_list, min_tree_depth=3):
        """
        Args:
            min_tree_depth (Optional[int]): minimum allowable depth for this taxonomy to be considered
        """
        filtered_list = []
        for taxonomy in taxonomy_list:
            tree_depth = 0
            try:
                current_taxonomy = self.tree[taxonomy].node_id
            except TypeError:
                # taxonomy represented in the reference database, but is not present in the tree
                continue
            while not current_taxonomy == "1":
                tree_depth += 1
                current_taxonomy = self.tree[current_taxonomy].parent_id
            if tree_depth < min_tree_depth:
                continue
            filtered_list.append(taxonomy)
        return filtered_list

    def taxonomic_lineage(self, taxonomy):
        """

        Args:
            taxonomy (str): taxonomy name or taxonomy ID

        Returns:
            list of lineage
        """
        taxonomy_id = self.tree[taxonomy].node_id
        lineage = [taxonomy_id]
        while not taxonomy_id == "1":
            taxonomy_id = self.tree[taxonomy_id].parent_id
            lineage.insert(0, taxonomy_id)
        return lineage

    def lca_majority(self, taxonomy_list, majority_cutoff):
        """Finds a consensus majority up a tree structure.

        Args:
            taxonomy_list (list): list of taxonomy names
            majority_cutoff (float): this is the total length of the taxonomy list * majority fraction

        Returns:
            list: string representing majority taxonomy ID, lineage counts dictionary with a key
                of taxonomy ID and value of list populated with ordered lineage taxonomy IDs

        Example:
            >>> taxonomy_list = ['gamma subgroup', 'RNA similarity group I',
                                 'purple photosynthetic bacteria and relatives',
                                 'not Bacteria Haeckel 1894',
                                 'purple photosynthetic bacteria and relatives', 'gamma subgroup',
                                 'gamma subgroup', 'purple photosynthetic bacteria and relatives',
                                 'purple photosynthetic bacteria and relatives']
            >>> majority_cutoff = len(taxonomy_list) * .5
            >>> tree = Tree("ref/ncbi_taxonomy_tree.txt")
            >>> tree.lca_majority(taxonomy_list)
            ('1224', {'1224': ['1', '131567', '99999999', '2', '1224'],
                      '1236': ['1', '131567', '99999999', '2', '1224', '1236'],
                      '2': ['1', '131567', '99999999', '2'],
                      '286': ['1', '131567', '99999999', '2', '1224', '1236', '72274', '135621', '286']})

        """
        lineage_counts = Counter()
        lineages = {}
        for taxonomy in taxonomy_list:
            lineage = self.taxonomic_lineage(taxonomy)
            lineage_counts.update(lineage)
            lineages[self.tree[taxonomy].node_id] = lineage
        lineage_indexes = index_of_list_items(lineages.values())
        for taxonomy in lineage_indexes:
            if lineage_counts[taxonomy] > majority_cutoff:
                return taxonomy, lineages
        return "1", lineages

    def counts_to_majority_list(self, taxonomy_counts, lineages, majority_id):
        aggregate_counts = []
        for taxonomy, taxonomy_count in taxonomy_counts.items():
            taxonomy_id = self.tree[taxonomy].node_id
            if majority_id in lineages[taxonomy_id]:
                taxonomy_id = majority_id
            aggregate_counts.extend([taxonomy_id] * taxonomy_count)
        return aggregate_counts

    def lca_star(self, taxonomy_list, min_tree_depth=3, majority_threshold=0.50):
        """

        Args:
            taxonomy_list (list):
            min_tree_depth (int):
            majority_threshold (float): 0-1;

        Returns:
            dict of 'taxonomy' and 'pvalue'
        """
        # tree depth based filter
        taxonomy_list = self.filter_taxonomy_list(taxonomy_list, min_tree_depth)
        # all have been filtered
        if not taxonomy_list:
            majority = "root"
            p = 1.
        else:
            taxonomy_counts = Counter(taxonomy_list)
            total_count = len(taxonomy_list)
            majority_cutoff = total_count * majority_threshold
            # majority based on existing taxonomy counts alone
            if taxonomy_counts.most_common()[0][1] > majority_cutoff:
                majority = taxonomy_counts.most_common()[0][0]
                p = nettleton_pvalue(taxonomy_list, majority)
            # create majority from lca
            else:
                majority, lineages = self.lca_majority(taxonomy_list, majority_cutoff)
                aggregate_counts = self.counts_to_majority_list(taxonomy_counts, lineages, majority)
                p = nettleton_pvalue(aggregate_counts, majority)
        return {"taxonomy":majority, "pvalue":p}


class Node(object):

    def __init__(self, taxonomy, node_id, parent_id):
        # the current node's string ID
        self.taxonomy = taxonomy
        # the current node's digit ID
        self.node_id = node_id
        self.parent_id = parent_id


class BlastHit(object):

    def __init__(self, taxonomies=None):
        """
        Args:
            taxonomies (Optional[list]): when initiated with a tax list; :func:`best_hit` and
                :func:`add` will no longer operate as intended
        """
        if taxonomies is None:
            # increasing bitscore sorted
            self.taxonomies = deque()
            self.bitscores = deque()
        else:
            self.taxonomies = taxonomies

    def __repr__(self):
        return "{cls}[{tax}]".format(cls=self.__class__.name, tax=self.taxonomies)

    def add(self, taxonomy, bitscore, max_hits=10, top_fraction=None):
        """Add entry to this :class:`BlastHit` group.

        Args:
            taxonomy (str): taxonomy name
            bitscore (str): bitscore for hit
            max_hits (int): maximum number of hits to consider for this :class:`BlastHit` group
            top_fraction (float): fraction cutoff from best bitscore, e.g. 0.3 will filter out 699 when best bitscore is 1000

        Notes:
            max_hits and top_fraction work in conjunction of one another
        """
        bitscore = float(bitscore)

        if top_fraction and self.bitscores:
             if bitscore < self.bitscores[-1] - (self.bitscores[-1] * top_fraction):
                 bitscore = None

        if bitscore:
            # insert into sorted list
            idx = bisect.bisect_left(self.bitscores, bitscore)
            self.bitscores.insert(idx, bitscore)
            self.taxonomies.insert(idx, taxonomy)
            if len(self.taxonomies) > max_hits:
                # remove lowest bitscore
                self.taxonomies.popleft()
                self.bitscores.popleft()

    def best_hit(self):
        return self.taxonomies[-1]

    def majority(self):
        # no repeated taxonomies
        if len(self.taxonomies) == len(set(self.taxonomies)):
            return self.best_hit()
        else:
            # count each taxonomy, grab top taxonomy
            most_common = Counter(self.taxonomies).most_common(1)[0][0]
            # need to flip to grab best bitscore
            self.taxonomies.reverse()
            # left most index match
            idx = self.taxonomies.index(most_common)
            return self.taxonomies[idx]


def index_of_list_items(lists):
    """Find the highest index for items among list of lists and return ordered dictionary sorted by
    decreasing index value.

    Args:
        lists (list): list of lists

    Returns:
        OrderedDict of taxonomy ID, sorted by depth with deepest node first

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


def parse_annotation_map(path):
    """Parse file and return dict of column 1 to value of column 2.

    Args:
        path (str): file path to annotation map file of name to product

    Returns:
        dict: name to product

    Notes:
        For RefSeq, the annotation map file would look like:

            >gi|494717612|ref|WP_007453478.1| hypothetical protein [Micromonospora lupini]
            >gi|494717623|ref|WP_007453489.1| Tryptophanyl-tRNA synthetase (fragment) [Micromonospora lupini]
            >gi|494717625|ref|WP_007453491.1| Tryptophanyl-tRNA synthetase (fragment) [Micromonospora lupini]

        Where the name and product are separated by a space.

    """
    logging.debug("Reading the annotation map")
    am = {}
    with open(path) as fh:
        for line in fh:
            toks = line.strip().partition(" ")
            am[toks[0].strip(">")] = toks[2] if toks[2] else "hypothetical protein"
    return am


def parse_ncbi_map(path):
    """Parse file and create dict of column 1 to value of column 2.

    Args:
        path (str): file path to NCBI map file of taxonomic ID to Taxonomic name

    Returns:
        dict: taxonomy ID number (as str) to taxonomy name

    Notes:
        The NCBI map file (tab delimited) looks like:

            1  root	         -1   0
            2  Bacteria	     -1   0
            6  Azorhizobium  -1  98

        Gives:

            dict('1':'root', '2':'Bacteria', '6':'Azorhizobium')

    """
    logging.debug("Parsing NCBI map file")
    m = {}
    with open(path) as fh:
        for line in fh:
            toks = line.strip().split("\t")
            m[toks[0]] = toks[1]
    return m


def parse_blast_results(blast_tab, annotation_map, orf_summary, tree=None,
                        min_identity=60, min_bitscore=0, min_length=60, max_evalue=0.000001,
                        max_hits_per_orf=25, top_fraction_of_hits=None, lca_threshold=1):
    """Parse BLAST results (-outfmt 6), filter, and aggregate ORF taxonomies.

    Args:
        blast_tab (str):
        annotation_map (str):
        orf_summary (dict):

        lca_threshold (float): the first parent above this fraction of representation (its count is
            greater than the total * lca_threshold)

    Returns:
        dict of dicts where first key is contig name, inner key is ORF ID with a value of its
            consensus taxonomy name as a string

    Raises:
        AssertionError when ORF summary method is not supported (['lca', 'best', 'majority'])
    """
    assert orf_summary in ['lca', 'best', 'majority']
    annotations = parse_annotation_map(annotation_map)

    if orf_summary == "lca" and tree is None:
        logging.critical("LCA summaries require a taxonomic tree")
        sys.exit(1)

    # ec_re = re.compile(r'(\d+[.]\d+[.]\d+[.]\d+)')
    tax_re = re.compile(r'\[([^\[]+)\]')
    blast_6 = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
               'sstart', 'send', 'evalue', 'bitscore']

    contigs = defaultdict(lambda: defaultdict(BlastHit))

    with open(blast_tab) as blast_tab_fh:
        current_hit = ""
        hit_count = 0

        for hsp in blast_tab_fh:
            toks = dict(zip(blast_6, hsp.strip().split("\t")))

            # user filtering cutoffs
            if int(toks['length']) < min_length or float(toks['pident']) < min_identity or \
               float(toks['bitscore']) < min_bitscore or float(toks['evalue']) > max_evalue:
                continue

            try:
                raw_product = annotations[toks['sseqid']]
            except KeyError:
                logging.critical("Annotations are missing from your map [see: %s] that appear in your BLAST reference database." % toks['sseqid'])
                sys.exit(1)
            # try:
            #     ec = ec_re.findall(raw_product)[0]
            # except IndexError:
            #     ec = ""
            try:
                taxonomy = tax_re.findall(raw_product)[0]
            except IndexError:
                taxonomy = ""

            contig_name, _, orf_idx = toks['qseqid'].rpartition("_")

            bbb = toks['bitscore']
            contigs[contig_name][orf_idx].add(taxonomy, toks['bitscore'], max_hits_per_orf, top_fraction_of_hits)

    # can't guarantee ordering, so iterate again and aggregate ORF assignments
    orf_assignments = defaultdict(dict)
    for contig, orfs in contigs.items():
        for orf, hits in orfs.items():
            if orf_summary == "best":
                orf_assignments[contig][orf] = hits.best_hit()
            elif orf_summary == "majority":
                orf_assignments[contig][orf] = hits.majority()
            # orf_summary == "lca":
            else:
                hits.taxonomies.reverse()
                orf_assignments[contig][orf] = tree.lca(hits.taxonomies, threshold=lca_threshold)
    return orf_assignments


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


def process_orfs(orf_assignments, tree, taxonomy_name_map=None):
    """Processing the already classified ORFs through secondary contig classification.

    Args:
        orf_assignments

    """
    print("contig", "lca_star", "lca_star_p", "majority", "majority_p", "lca_squared", sep="\t")
    for contig, orfs in orf_assignments.items():
        taxonomies = list(orfs.values())
        lca_star_result = tree.lca_star(taxonomies)
        second_lca = tree.lca(taxonomies, threshold=1)
        majority = BlastHit(taxonomies).majority()
        majority_p = nettleton_pvalue(taxonomies, majority)

        if taxonomy_name_map:
            lca_star_lineage = ";".join([tree.taxonomy_id_to_name(taxonomy_name_map, i) for i in tree.taxonomic_lineage(lca_star_result['taxonomy'])])
            second_lca_lineage = ";".join([tree.taxonomy_id_to_name(taxonomy_name_map, i) for i in tree.taxonomic_lineage(second_lca)])
            majority_lineage = ";".join([tree.taxonomy_id_to_name(taxonomy_name_map, i) for i in tree.taxonomic_lineage(majority)])
        else:
            lca_star_lineage = ";".join(["%s (%s)" % (tree[i].taxonomy, i) for i in tree.taxonomic_lineage(lca_star_result['taxonomy'])])
            second_lca_lineage = ";".join(["%s (%s)" % (tree[i].taxonomy, i) for i in tree.taxonomic_lineage(second_lca)])
            majority_lineage = ";".join(["%s (%s)" % (tree[i].taxonomy, i) for i in tree.taxonomic_lineage(majority)])
        print(contig, lca_star_lineage, lca_star_result['pvalue'], majority_lineage, majority_p, second_lca_lineage, sep="\t")


# def main(fasta, blasttab, annotation_map, out_file, blast_lambda, blast_k, min_bsr, min_identity, min_bitscore, min_length, max_evalue, max_hits):
    # logging.info("Finding best possible alignment score per sequence")
annotation_map = "/Users/brow015/devel/blastlca/refseq-nr-2014-01-18-names.txt"
taxonomy_tree = "/Users/brow015/devel/blastlca/ncbi_taxonomy_tree.txt"
tree = Tree(taxonomy_tree)
# best_scores = calculate_best_score(fasta, blast_lambda, blast_k)
blasttab = "/Users/brow015/devel/blastlca/Metat_1000000_contigs.refseq-nr-2014-01-18.diamond"
taxonomy_name_map = parse_ncbi_map("/Users/brow015/devel/blastlca/ncbi.map")
orf_assignments = parse_blast_results(blasttab, annotation_map, orf_summary="lca", tree=tree, min_identity=60, min_bitscore=0, min_length=60, max_evalue=0.000001, max_hits_per_orf=10, top_fraction_of_hits=1)
process_orfs(orf_assignments, tree, taxonomy_name_map=taxonomy_name_map)
    # logging.info("Complete")
    # out_file = "test.parsed.new"
    # ncbi_megan_map = "ncbi.map"
    # tree_file = "ncbi_taxonomy_tree.txt"
    # orf_summary = "majority"

    # ncbi_tree = LCAStar(ncbi_tree, min_depth=1, alpha=0.51, min_reads=1)
    ## calculate taxonomy statistics LCA,  for each ORF
    # contig_to_taxa = {}
    ## LCA^2, Majority, and LCA* for each ORF
    # writeout(args, contig_to_lca, lcastar, ncbi_megan_map)


# if __name__ == "__main__":
#
#     def _file_exists(parser, arg):
#         if not os.path.exists(arg):
#             parser.error("The file %s does not exist" % arg)
#         if not os.path.isfile(arg):
#             parser.error("Expected file, not folder (%s)" % arg)
#         return arg
#
#     p = argparse.ArgumentParser(description=__doc__,
#             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#     p.add_argument("fasta", type=lambda x: _file_exists(p, x), help="query FASTA file path")
#     p.add_argument("blasttab", type=lambda x: _file_exists(p, x), help="BLAST (-outfmt 6) result file path")
#     p.add_argument("map", type=lambda x: _file_exists(p, x), help="BLAST database name map with space separated name product per line")
#     p.add_argument("outfile", help="parsed blast results")
#     # p.add_argument("--blast-lambda", type=float, default=0.267, help="lambda parameter for custom matrix")
#     # p.add_argument("--blast-k", type=float, default=0.041, help="k parameter for custom matrix")
#     # p.add_argument("--min-bsr", type=float, default=0.4, help="minimum allowable bitscore ratio (the bitscore of a query and target sequence over the bitscore when both query and target are the query sequence)")
#     p.add_argument("--min-identity", type=float, default=60, help="minimum allowable percent ID of hsp")
#     p.add_argument("--min-bitscore", type=float, default=0, help="minimum allowable bitscore of hsp")
#     p.add_argument("--min-length", type=int, default=60, help="minimum allowable alignment length of hsp")
#     p.add_argument("--max-evalue", type=float, default=0.000001, help="maximum allowable evalue of hsp")
#     p.add_argument("--max-hits", type=int, default=25, help="maximum hits per query to be considered for LCA or majority")
#     args = p.parse_args()
#     logging.basicConfig(level=logging.INFO, datefmt="%Y-%m-%d %H:%M", format="[%(asctime)s] %(message)s")
#     main(args.fasta, args.blasttab, args.map, args.outfile, args.blast_lambda, args.blast_k, args.min_bsr, args.min_identity, args.min_bitscore, args.min_length, args.max_evalue, args.max_hits)
