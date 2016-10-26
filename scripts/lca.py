#!/usr/bin/env python
# coding=utf-8

import bisect
import click
import logging
import math
import os
import re
import sys
from collections import Counter, defaultdict, deque, OrderedDict
from math import log, sqrt, erfc


class Node(object):

    def __init__(self, taxonomy, node_id, parent_id):
        """Represents a node within a tree.

        Args:
            taxonomy (str): taxonomy name or ID
            node_id (str): taxonomy ID
            parent_id (str): taxonomy ID of parent

        """
        # the current node's string ID
        self.taxonomy = taxonomy
        # the current node's digit ID
        self.node_id = node_id
        self.parent_id = parent_id


class Tree(object):

    def __init__(self, tax_tree):
        """Builds reference dictionary from tab delimited text of Taxonomy Name, Taxonomy ID,
        Parent Taxonomy ID.

        Args:
            tax_tree (str): file path to taxonomy tree txt

        Raises:
            AssertionError when child and parent IDs are equal when not root

        Notes:
            An example of the file:

                root	    1	        1
                all	        1	        1
                prokaryotes	99999999    131567

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
        """Adds node to tree dictionary.

        Args:
            taxonomy (str): the taxonomy name
            node_id (str): the taxonomy id
            parent_id (str): the parent's taxonomy id

        """
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
            taxonomies (list): list of taxonomies (ID or name); when using threshold < 1 they
                should be ordered by decreasing bitscore
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
        """Filters a taxonomy list by tree depth in an effort to classify at a higher resolution.

        Args:
            taxonomy_list (list): list of taxonomy names or IDs to filter
            min_tree_depth (Optional[int]): minimum allowable depth for this taxonomy to be considered

        Returns:
            list

        Example:
            >>> tree = Tree("ref/ncbi_taxonomy_tree.txt")
            >>> tree.filter_taxonomy_list(["bacteria", "Pseudomonadaceae", "Prokaryotae"], min_tree_depth=1)
            ['bacteria', 'Pseudomonadaceae', 'Prokaryotae']
            >>> tree.filter_taxonomy_list(["bacteria", "Pseudomonadaceae", "Prokaryotae"], min_tree_depth=3)
            ['bacteria', 'Pseudomonadaceae']
            >>> tree.filter_taxonomy_list(["bacteria", "Pseudomonadaceae", "Prokaryotae"], min_tree_depth=4)
            ['Pseudomonadaceae']

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
        """For a given taxonomy name or ID, return its lineage as a list of IDs.

        Args:
            taxonomy (str): taxonomy name or taxonomy ID

        Returns:
            list of lineage

        Example:
            >>> tree = Tree("ref/ncbi_taxonomy_tree.txt")
            >>> tree.taxonomic_lineage("Pseudomonadaceae")
            ['1', '131567', '99999999', '2', '1224', '1236', '72274', '135621']

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
            taxonomy_list (list): list of taxonomy names or IDs
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
        """Aggregate the counts across lineage observations for the majority ID.

        Args:
            taxonomy_counts (collections.Counter): count per taxon
            lineages (list): list of lineages per taxon
            majority_id (str): the taxonomy name or ID upon which to aggregate counts

        Returns:
            list of representative taxonomies

        Example:
            >>> tree = Tree("ref/ncbi_taxonomy_tree.txt")
            >>> taxonomy_list = ['gamma subgroup', 'RNA similarity group I',
                                 'purple photosynthetic bacteria and relatives',
                                 'not Bacteria Haeckel 1894',
                                 'purple photosynthetic bacteria and relatives', 'gamma subgroup',
                                 'gamma subgroup', 'purple photosynthetic bacteria and relatives',
                                 'purple photosynthetic bacteria and relatives']
            >>> majority_id, lineages = tree.lca_majority(taxonomy_list, 0.5 * len(taxonomy_list))
            >>> tree.counts_to_majority_list(Counter(taxonomy_list), lineages, majority_id)
            ['1224', '2', '1224', '1224', '1224', '1224', '1224', '1224', '1224']

        """
        aggregate_counts = []
        for taxonomy, taxonomy_count in taxonomy_counts.items():
            taxonomy_id = self.tree[taxonomy].node_id
            if majority_id in lineages[taxonomy_id]:
                taxonomy_id = majority_id
            aggregate_counts.extend([taxonomy_id] * taxonomy_count)
        return aggregate_counts

    def lca_star(self, taxonomy_list, min_tree_depth=3, majority_threshold=0.50):
        """Find the LCA within a list of taxonomies after filtering those taxonomies by tree depth.
        One can also vary what constitutes a majority consensus for the counts, with the default
        being 50%.

        Args:
            taxonomy_list (list): list of taxonomy names or IDs
            min_tree_depth (int): the mininum allowable tree depth of taxon to be considered within
                the taxonomy list; those found sooner in the tree will be filtered out of consideration
            majority_threshold (float): 0-1; the fraction of taxonomy counts which constitutes a
                majority; a lower fraction will classify with less confidence deeper in the tree
                while a higher threshold will classify with more confidence higher in the tree

        Returns:
            dict of 'taxonomy' and 'pvalue'

        Example:
            >>> tree = Tree("ref/ncbi_taxonomy_tree.txt")
            >>> taxonomy_list = ['gamma subgroup', 'RNA similarity group I',
                                 'purple photosynthetic bacteria and relatives',
                                 'not Bacteria Haeckel 1894',
                                 'purple photosynthetic bacteria and relatives', 'gamma subgroup',
                                 'gamma subgroup', 'purple photosynthetic bacteria and relatives',
                                 'purple photosynthetic bacteria and relatives']
            >>> tree.lca_star(taxonomy_list)
            {'pvalue': 0.012791848981090311, 'taxonomy': '1224'}

        """
        # tree depth based filter
        taxonomy_list = self.filter_taxonomy_list(taxonomy_list, min_tree_depth)
        # all have been filtered
        if not taxonomy_list:
            majority = "root"
            p = 1.
        else:
            taxonomy_counts = Counter(taxonomy_list)
            majority_cutoff = len(taxonomy_list) * majority_threshold
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


class BlastHits(object):

    def __init__(self, names=None, max_hits=10, top_fraction=None):
        """Class that represents BLAST hits for a single target sequence. Hits are added to queues
        for bitscore and ID and ordered by increasing bitscore.

        Args:
            names (Optional[list]): when initiated with a name list; :func:`best_hit` and
                :func:`add` will no longer operate as intended
            max_hits (int): maximum number of hits to consider for this :class:`BlastHits` group
            top_fraction (float): fraction cutoff from best bitscore, e.g. 0.3 will filter out 699 when best bitscore is 1000

        Notes:
            max_hits and top_fraction work in conjunction of one another

        Examples:
            >>> hits = BlastHits()
            >>> hits.add("name_1", 100)
            >>> hits.add("name_2", 110)
            >>> hits.add("name_3", 120)
            >>> hits.add("name_4", 130)
            >>> hits.add("name_4", 140)
            >>> hits.majority()
            'name_4'
            >>> hits.best_hit()
            'name_4'
            >>> hits.add("name_5", 150)
            >>> hits.best_hit()
            'name_5'
            # demo max hits and top fraction
            >>> hits = BlastHits(max_hits=3, top_fraction=0.7)
            >>> hits.add("name_1", 100)
            # filtered hit
            >>> hits.add("name_2", 65)
            >>> hits.names
            deque(['name_1'])
            >>> hits.add("name_3", 70)
            deque(['name_3', 'name_1'])
            >>> hits.add("name_4", 80)
            >>> hits.add("name_5", 85)
            >>> hits.names
            deque(['name_4', 'name_5', 'name_1'])
            # new best bitscore
            >>> hits.add("name_6", 125)
            >>> dict(zip(hits.names, hits.bitscores))
            {'name_1': 100.0, 'name_6': 125.0}

        """
        if names is None:
            # increasing bitscore sorted
            self.names = deque()
            self.bitscores = deque()
        else:
            self.names = names
        self.max_hits = max_hits
        self.top_fraction = top_fraction

    def __repr__(self):
        return "{cls}[{tax}]".format(cls=self.__class__.__name__, tax=self.names)

    def add(self, name, bitscore):
        """Add entry to this :class:`BlastHits` group.

        Args:
            name (str): hit identifier
            bitscore (str): bitscore for hit

        """
        bitscore = float(bitscore)

        if self.top_fraction and self.bitscores:
            # the filter
            if bitscore < (self.bitscores[-1] * self.top_fraction):
                bitscore = None
            # new best
            elif bitscore > self.bitscores[-1]:
                score = self.bitscores[0]
                while score < bitscore * self.top_fraction:
                    self.names.popleft()
                    self.bitscores.popleft()
                    score = self.bitscores[0]

        if bitscore:
            # insert into sorted list
            idx = bisect.bisect_left(self.bitscores, bitscore)
            self.bitscores.insert(idx, bitscore)
            self.names.insert(idx, name)
            if len(self.names) > self.max_hits:
                # remove lowest bitscore
                self.names.popleft()
                self.bitscores.popleft()

    def best_hit(self):
        """Returns the hit ID of the best scoring alignment."""
        return self.names[-1]

    def majority(self):
        """Returns the hit ID of the best scoring hit ID that is repeated or the best hit when
        no items are repeated.
        """
        # no repeated names
        if len(self.names) == len(set(self.names)):
            return self.best_hit()
        else:
            # count each taxonomy, grab top taxonomy
            most_common = Counter(self.names).most_common(1)[0][0]
            # need to flip to grab best bitscore
            names_reversed = self.names.copy()
            names_reversed.reverse()
            # left most index match
            idx = names_reversed.index(most_common)
            return names_reversed[idx]


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option("0.0.0")
@click.pass_context
def cli(obj):
    "LCA methods for BLAST tabular output."


@cli.command("tree-based", short_help="enables tree based LCA and LCA star methods")
@click.argument("tsv", type=click.Path(exists=True))
@click.argument("annotation", type=click.Path(exists=True))
@click.argument("tree", type=click.Path(exists=True))
@click.argument("output", type=click.File("w"))
@click.option("-m", "--name-map", type=click.Path(exists=True), help="taxonomy name map that uses taxonomy IDs within TREE and converts them to the name in this file, e.g. 'ncbi.map'; taxonomy IDs in column 1, new taxonomy name in column 2")
@click.option("-s", "--summary-method", type=click.Choice(["lca", "majority", "best"]), default="lca", show_default=True, help="summary method for annotating ORFs")
@click.option("--min-identity", type=int, default=60, show_default=True, help="minimum allowable percent ID of BLAST hit")
@click.option("--min-bitscore", type=int, default=0, show_default=True, help="minimum allowable bitscore of BLAST hit; 0 disables")
@click.option("--min-length", type=int, default=60, show_default=True, help="minimum allowable BLAST alignment length")
@click.option("--max-evalue", type=float, default=0.000001, show_default=True, help="maximum allowable e-value of BLAST hit")
@click.option("--max-hits", type=int, default=10, show_default=True, help="maximum number of BLAST hits to consider when summarizing ORFs; can drastically alter LCA assignments if too high; when too low, will affect error probability calculation resulting in higher p-values")
@click.option("--top-fraction", type=float, default=1, show_default=True, help="filters ORF BLAST hits by only keep hits within this fraction of the highest bitscore")
def tree_based_parsing(tsv, annotation, tree, output, name_map, summary_method, min_identity,
                       min_bitscore, min_length, max_evalue, max_hits, top_fraction):
    """Parse TSV (tabular BLAST output [-outfmt 6]), grabbing taxonomy from ANNOTATION, and use
    TREE to compute LCAs.

    Annotation file should include your BLAST subject sequence ID, a space, then product:

        \b
        >id1 hypothetical protein [Micromonospora lupini]
        >id2 conserved exported hypothetical protein [Micromonospora lupini]

    Taxonomy is contained with brackets '[]'.

    The optional name map, takes derived taxonomy IDs from TREE and renames the taxonomy to its
    match within this file. An example header for this tab-delimited file looks like:

        \b
        9	Buchnera aphidicola
        10	Cellvibrio
        11	[Cellvibrio] gilvus
        13	Dictyoglomus

    """
    logging.basicConfig(level=logging.INFO, datefmt="%Y-%m-%d %H:%M", format="[%(asctime)s] %(message)s")
    logging.info("Parsing %s" % tree)
    tree = Tree(tree)
    logging.info("Parsing %s" % tsv)
    orf_assignments = parse_blast_results(tsv, annotation, orf_summary=summary_method, tree=tree,
                                          min_identity=min_identity, min_bitscore=min_bitscore,
                                          min_length=min_length, max_evalue=max_evalue,
                                          max_hits_per_orf=max_hits, top_fraction_of_hits=top_fraction)
    logging.info("Assigning taxonomies to contigs")
    # TODO: this should include name mapped 'product' column so we propogate function
    process_orfs(orf_assignments, tree, output, taxonomy_name_map=name_map)
    logging.info("Complete")


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
    logging.debug("Reading the annotation map [%s]" % path)
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
    logging.debug("Parsing NCBI map file [%s]" % path)
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

    contigs = defaultdict(lambda: defaultdict(lambda: BlastHits(max_hits=max_hits_per_orf,
                                                                top_fraction=top_fraction_of_hits)))

    with open(blast_tab) as blast_tab_fh:
        current_hit = ""
        hit_count = 0

        for hsp in blast_tab_fh:
            toks = dict(zip(blast_6, hsp.strip().split("\t")))

            # user filtering cutoffs
            if int(toks['length']) < min_length or float(toks['pident']) < min_identity or float(toks['evalue']) > max_evalue: continue
            if min_bitscore and float(toks['bitscore']) < min_bitscore: continue

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
                # this makes this a non-generic blast parser
                taxonomy = tax_re.findall(raw_product)[0]
                contig_name, _, orf_idx = toks['qseqid'].rpartition("_")
                contigs[contig_name][orf_idx].add(taxonomy, toks['bitscore'])
            except IndexError:
                continue

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
                hits.names.reverse()
                orf_assignments[contig][orf] = tree.lca(hits.names, threshold=lca_threshold)
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


def process_orfs(orf_assignments, tree, output, taxonomy_name_map=None):
    """Processing the already classified ORFs through secondary contig classification.

    Args:
        orf_assignments (dict): dict of dict for per ORF tax assignment per contig
        tree (Tree): taxonomic tree object
        output (filehandle): output file handle
        taxonomy_name_map (str): file path to taxonomy name map
    """
    if taxonomy_name_map:
        taxonomy_name_map = parse_ncbi_map(taxonomy_name_map)

    print("contig", "lca_star", "lca_star_p", "majority", "majority_p", "lca_squared", sep="\t",
          file=output)
    for contig, orfs in orf_assignments.items():
        taxonomies = list(orfs.values())
        lca_star_result = tree.lca_star(taxonomies)
        second_lca = tree.lca(taxonomies, threshold=1)
        majority = BlastHits(taxonomies).majority()
        majority_p = nettleton_pvalue(taxonomies, majority)

        if taxonomy_name_map:
            lca_star_lineage = ";".join([tree.taxonomy_id_to_name(taxonomy_name_map, i) for i in tree.taxonomic_lineage(lca_star_result['taxonomy'])])
            second_lca_lineage = ";".join([tree.taxonomy_id_to_name(taxonomy_name_map, i) for i in tree.taxonomic_lineage(second_lca)])
            majority_lineage = ";".join([tree.taxonomy_id_to_name(taxonomy_name_map, i) for i in tree.taxonomic_lineage(majority)])
        else:
            lca_star_lineage = ";".join(["%s (%s)" % (tree[i].taxonomy, i) for i in tree.taxonomic_lineage(lca_star_result['taxonomy'])])
            second_lca_lineage = ";".join(["%s (%s)" % (tree[i].taxonomy, i) for i in tree.taxonomic_lineage(second_lca)])
            majority_lineage = ";".join(["%s (%s)" % (tree[i].taxonomy, i) for i in tree.taxonomic_lineage(majority)])
        print(contig, lca_star_lineage, lca_star_result['pvalue'], majority_lineage, majority_p,
              second_lca_lineage, sep="\t", file=output)


if __name__ == "__main__":
    cli()
