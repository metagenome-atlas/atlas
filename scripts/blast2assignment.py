#!/usr/bin/env python
# coding=utf-8
import bisect
import click
import contextlib
import gzip
import logging
import math
import os
import re
import shelve
import sys
from collections import Counter, defaultdict, deque, OrderedDict
from itertools import groupby
from math import log, sqrt, erfc


BLAST6 = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
          "qend", "sstart", "send", "evalue", "bitscore"]


logging.basicConfig(level=logging.INFO, datefmt="%Y-%m-%d %H:%M", format="[%(asctime)s] %(message)s")
gzopen = lambda f: gzip.open(f, mode="rt") if f.endswith(".gz") else open(f)


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

    def __init__(self, tree_file):
        """Builds reference dictionary of Taxonomy Name, Taxonomy ID, and Parent Taxonomy ID."""
        self.tree = defaultdict(dict)

        with open(tree_file) as tree_fh:
            for line in tree_fh:
                toks = line.strip().split("\t")
                if not toks[0] == '1' and not toks[2] == '1':
                    assert not toks[0] == toks[2]
                if not len(toks) == 3:
                    logging.warning("Line [%s] does not have ID, NAME, PARENTID" % line.strip())
                    continue
                self.add_node(toks[1], toks[0], toks[2])

    def add_node(self, taxonomy, node_id, parent_id):
        """Adds node to tree dictionary.

        Args:
            taxonomy (str): the taxonomy name
            node_id (str): the taxonomy id
            parent_id (str): the parent's taxonomy id

        """
        # taxonomy id to node mapping; ensures unique nodes despite non-unique names
        self.tree[node_id] = Node(taxonomy, node_id, parent_id)
        # taxonomy string to node mapping
        # self.tree[taxonomy] = self.tree[node_id]

    def taxonomy_id_to_name(self, key):
        """
        Args:
            key (str): taxonomy ID used to search `id_map`

        Returns:
            string of mapped key along with key, e.g. New Name (1224)
        """
        return "{name} ({digit})".format(name=self.tree[key].taxonomy, digit=key)

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
                    return self.tree[current_taxonomy].node_id
                # traverse up tree
                current_taxonomy = self.tree[current_taxonomy].parent_id
        return "1"

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
            current_taxonomy = taxonomy
            # try:
            #     current_taxonomy = self.tree[taxonomy].node_id
            # except TypeError:
            #     # taxonomy represented in the reference database, but is not present in the tree
            #     continue
            while not current_taxonomy == "1":
                tree_depth += 1
                current_taxonomy = self.tree[current_taxonomy].parent_id
            if tree_depth < min_tree_depth:
                continue
            filtered_list.append(taxonomy)
        return filtered_list

    def taxonomic_lineage(self, taxonomy):
        """For a given taxonomy ID, return its lineage as a list of IDs.

        Args:
            taxonomy (str): taxonomy name or taxonomy ID

        Returns:
            list of lineage

        Example:
            >>> tree = Tree("ref/ncbi_taxonomy_tree.txt")
            >>> tree.taxonomic_lineage("135621")
            ['1', '131567', '99999999', '2', '1224', '1236', '72274', '135621']

        """
        # a large portion of ORFs
        if taxonomy == "1":
            return [taxonomy]

        lineage = [taxonomy]
        while not taxonomy == "1":
            taxonomy = self.tree[taxonomy].parent_id
            # prepend
            lineage.insert(0, taxonomy)
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
            lineages[taxonomy] = lineage
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
            taxonomy_id = taxonomy
            if majority_id in lineages[taxonomy_id]:
                taxonomy_id = majority_id
            aggregate_counts.extend([taxonomy_id] * taxonomy_count)
        return aggregate_counts

    def lca_star(self, taxonomy_list, min_tree_depth=3, majority_threshold=0.51):
        """Find the LCA within a list of taxonomies after filtering those taxonomies by tree depth.
        One can also vary what constitutes a majority consensus for the counts, with the default
        being 51%.

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
            majority = "1"
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
            top_fraction (float): fraction cutoff from best bitscore, e.g. 0.7 will filter out 699 when best bitscore is 1000

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

    def __len__(self):
        return len(self.names)

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
@click.version_option("0.2.0")
@click.pass_context
def cli(obj):
    "Methods for BLAST tabular output."


@cli.command("eggnog", short_help="process blast hits and aggregate across ORFs")
@click.argument("tsv", type=click.Path(exists=True))
@click.argument("namemap", type=click.Path(exists=True))
@click.argument("output", type=click.File("w"))
@click.option("--summary-method", type=click.Choice(["majority", "best"]), default="best", show_default=True, help="summary method for annotating ORFs; when majority and there is no majority, best is used")
@click.option("--min-identity", type=int, default=60, show_default=True, help="minimum allowable percent ID of BLAST hit")
@click.option("--min-bitscore", type=int, default=0, show_default=True, help="minimum allowable bitscore of BLAST hit; 0 disables")
@click.option("--min-length", type=int, default=60, show_default=True, help="minimum allowable BLAST alignment length")
@click.option("--max-evalue", type=float, default=0.000001, show_default=True, help="maximum allowable e-value of BLAST hit")
@click.option("--top-fraction", type=float, default=1, show_default=True, help="filters ORF BLAST hits before finding majority by only keep hits within this fraction of the highest bitscore; this is recommended over --max-hits")
@click.option("--max-hits", type=int, default=10, show_default=True, help="maximum number of BLAST hits to consider when summarizing ORFs as a majority")
@click.option("--table-name", default="eggnog", help="table name within namemap database; expected columns are 'eggnog_ssid_b', 'uniprot_id', 'ko_id', 'kegg_id', 'kegg_ec', 'cog_id', 'cog_func_id', 'cog_product', 'cog_path'")
def eggnog_parsing(tsv, namemap, output, summary_method, min_identity, min_bitscore, min_length,
                   max_evalue, top_fraction, max_hits, table_name):
    """Parse BLAST hits from EGGNOG.

    The BLAST hits are assumed to be sorted by query with decreasing bitscores (best alignment first):

        \b
        sort -k1,1 -k12,12rn tsv > sorted_tsv

    """
    import sqlite3
    logging.info("Parsing %s" % tsv)

    print("contig", "orf", "uniprot_id", "ko_id", "kegg_id", "kegg_ec", "cog_id", "cog_func_id",
          "cog_product", "cog_path", "%s_evalue" % table_name, "%s_bitscore" % table_name,
          sep="\t", file=output)
    with contextlib.closing(sqlite3.connect(namemap)) as conn, gzopen(tsv) as blast_tab_fh:
        cursor = conn.cursor()
        for query, qgroup in groupby(blast_tab_fh, key=lambda x: x.partition("\t")[0]):

            contig_name, _, orf_idx = query.rpartition("_")
            hit_id = ""
            bitscore = "NA"
            evalue = "NA"
            uniprot_id = "NA"
            ko_id = "NA"
            kegg_id = "NA"
            kegg_ec = "NA"
            cog_id = "NA"
            cog_func_id = "NA"
            cog_product = "NA"
            cog_path = "NA"
            orf_hits = BlastHits(max_hits=max_hits, top_fraction=top_fraction)
            lines = []
            idx = 0

            # iterate over blast hits per ORF
            for hsp in qgroup:
                toks = dict(zip(BLAST6, hsp.strip().split("\t")))
                # legacy for files mapped with incorrect reference
                toks["sseqid"] = toks["sseqid"].partition(".")[-1]
                if (int(toks["length"]) < min_length or
                        float(toks["pident"]) < min_identity or
                        float(toks["evalue"]) > max_evalue):
                    continue
                if min_bitscore and float(toks["bitscore"]) < min_bitscore:
                    # input is sorted by decreasing bitscore
                    break

                if summary_method == "best":
                    hit_id = toks["sseqid"]
                    bitscore = toks["bitscore"]
                    evalue = toks["evalue"]
                    break

                orf_hits.add(toks["sseqid"] + "_%d" % idx, toks["bitscore"])
                idx += 1
                lines.append(toks)

            # summary method is majority and we have passing HSPs
            if not hit_id and lines:
                majority_id = orf_hits.majority()
                line_idx = int(majority_id.rpartition("_")[-1])
                toks = lines[line_idx]
                hit_id = toks["sseqid"]
                bitscore = toks["bitscore"]
                evalue = toks["evalue"]

            if hit_id:
                cursor.execute('SELECT uniprot_id, ko_id, kegg_id, kegg_ec, cog_id, cog_func_id, cog_product, cog_path FROM %s WHERE eggnog_ssid_b="%s"' % (table_name, hit_id))
                try:
                    uniprot_id, ko_id, kegg_id, kegg_ec, cog_id, cog_func_id, cog_product, cog_path = cursor.fetchone()
                # legacy before database was pruned; can have hits not in metadata
                except TypeError:
                    pass

            # print for this query
            print(contig_name, "%s_%s" % (contig_name, orf_idx), uniprot_id, ko_id, kegg_id,
                  kegg_ec, cog_id, cog_func_id, cog_product, cog_path, evalue, bitscore, sep="\t",
                  file=output)
    logging.info("Complete")


@cli.command("refseq", short_help="enables tree based LCA and LCA star methods")
@click.argument("tsv", type=click.Path(exists=True))
@click.argument("namemap", type=click.Path(exists=True))
@click.argument("treefile", type=click.Path(exists=True))
@click.argument("output", type=click.File("w"))
@click.option("-s", "--summary-method", type=click.Choice(["lca", "majority", "best"]), default="lca", show_default=True, help="summary method for annotating ORFs; when using LCA, it's recommended that one limits the number of hits using --top-fraction; 'best' is fastest")
@click.option("-a", "--aggregation-method", type=click.Choice(["lca", "lca-majority", "majority"]), default="lca-majority", show_default=True, help="summary method for aggregating ORF taxonomic assignments to contig level assignment; 'lca' will result in most stringent, least specific assignments")
@click.option("--majority-threshold", type=float, default=0.51, show_default=True, help="constitutes a majority fraction at tree node for 'lca-majority' ORF aggregation method")
@click.option("--min-identity", type=int, default=60, show_default=True, help="minimum allowable percent ID of BLAST hit")
@click.option("--min-bitscore", type=int, default=0, show_default=True, help="minimum allowable bitscore of BLAST hit; 0 disables")
@click.option("--min-length", type=int, default=60, show_default=True, help="minimum allowable BLAST alignment length")
@click.option("--max-evalue", type=float, default=0.000001, show_default=True, help="maximum allowable e-value of BLAST hit")
@click.option("--max-hits", type=int, default=10, show_default=True, help="maximum number of BLAST hits to consider when summarizing ORFs; can drastically alter ORF LCA assignments if too high without further limits")
@click.option("--table-name", default="refseq", help="table name within namemap database; expected columns are 'name', 'function', and 'taxonomy'")
@click.option("--top-fraction", type=float, default=1, show_default=True, help="filters ORF BLAST hits by only keep hits within this fraction of the highest bitscore; this is recommended over --max-hits")
def refseq_parsing(tsv, namemap, treefile, output, summary_method, aggregation_method,
                   majority_threshold, min_identity, min_bitscore, min_length, max_evalue,
                   max_hits, table_name, top_fraction):
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

    """
    logging.info("Parsing %s" % tsv)
    tree = Tree(treefile)
    orf_assignments = parse_blast_results_with_tree(tsv, namemap, orf_summary=summary_method, tree=tree,
                          min_identity=min_identity, min_bitscore=min_bitscore,
                          min_length=min_length, max_evalue=max_evalue,
                          max_hits_per_orf=max_hits, table_name=table_name,
                          top_fraction_of_hits=top_fraction)
    logging.info("Assigning taxonomies to contigs using %s" % aggregation_method)
    process_orfs_with_tree(orf_assignments, tree, output, aggregation_method, majority_threshold, table_name)
    logging.info("Complete")


def parse_blast_results_with_tree(blast_tab, name_map, orf_summary, tree, min_identity=60,
                                  min_bitscore=0, min_length=60, max_evalue=0.000001,
                                  max_hits_per_orf=10, top_fraction_of_hits=None,
                                  table_name="refseq", lca_threshold=1):
    """Parse BLAST results (-outfmt 6), filter, and aggregate ORF taxonomies.

    Args:
        blast_tab (str): file path to blast TSV file
        name_map (dict): dict of tuples from parse_tree_annotation
        orf_summary (dict): method of ORF annotation selection

        lca_threshold (float): the first parent above this fraction of representation (its count is
            greater than the total * lca_threshold)

    Returns:
        dict: dict of dicts where first key is contig name, inner key is ORF ID; values are tuple of
              protein function, taxonomy ID, bitscore, evalue

    Raises:
        AssertionError when ORF summary method is not supported (['lca', 'best', 'majority'])
    """
    import sqlite3

    assert orf_summary in ["lca", "best", "majority"]
    contigs = defaultdict(dict)

    with contextlib.closing(sqlite3.connect(name_map)) as conn, gzopen(blast_tab) as blast_tab_fh:
        cursor = conn.cursor()
        for query, qgroup in groupby(blast_tab_fh, key=lambda x: x.partition("\t")[0]):

            contig_name, _, orf_idx = query.rpartition("_")
            protein_function = "hypothetical protein"
            protein_set = False
            taxonomy_id = "1"
            bitscore = "NA"
            evalue = "NA"
            orf_hits = BlastHits(max_hits=max_hits_per_orf, top_fraction=top_fraction_of_hits)

            # iterate over blast hits per ORF
            for i, hsp in enumerate(qgroup):
                toks = dict(zip(BLAST6, hsp.strip().split("\t")))
                if (int(toks["length"]) < min_length or
                        float(toks["pident"]) < min_identity or
                        float(toks["evalue"]) > max_evalue):
                    continue
                if min_bitscore and float(toks["bitscore"]) < min_bitscore:
                    # input is sorted by decreasing bitscore
                    break

                cursor.execute('SELECT function, taxonomy FROM %s WHERE name="%s"' % (table_name, toks["sseqid"]))
                current_function, current_taxonomy = cursor.fetchone()

                # update taxonomy based on pident
                # current_taxonomy = tree.climb_tree(current_taxonomy, float(toks["pident"]))

                # protein function is always assigned to best, passing alignment
                if not protein_set:
                    protein_function = current_function
                    bitscore = toks["bitscore"]
                    evalue = toks["evalue"]
                    protein_set = True
                    if orf_summary == "best":
                        taxonomy_id = current_taxonomy
                        break
                orf_hits.add(current_taxonomy, toks["bitscore"])

                # TODO if majority, function should be grabbed from its sseqid

            # ensure we have passing hits
            if len(orf_hits) > 0 and not orf_summary == "best":
                if orf_summary == "majority":
                    taxonomy_id = orf_hits.majority()
                # perform LCA to obtain taxonomy ID
                else:
                    orf_hits.names.reverse()
                    taxonomy_id = tree.lca(orf_hits.names, threshold=lca_threshold)
            contigs[contig_name][orf_idx] = (protein_function, taxonomy_id, bitscore, evalue)

    return contigs


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


def process_orfs_with_tree(orf_assignments, tree, output, aggregation_method, majority_threshold=0.51, table_name="refseq"):
    """Processing the already classified ORFs through secondary contig classification.

    Args:
        orf_assignments (dict): dict of dict for per ORF tax assignment per contig
        tree (Tree): taxonomic tree object
        output (filehandle): output file handle
        aggregation_method (str): lca, lca-majority, or majority
        majority_threshold (float): constitutes a majority fraction at tree node for 'lca-majority' ORF aggregation method
    """
    print("contig", "orf", "taxonomy", "erfc", "orf_taxonomy", "%s_product" % table_name,
          "%s_evalue" % table_name, "%s_bitscore" % table_name, sep="\t", file=output)
    for contig, orfs in orf_assignments.items():
        taxonomies = [x[1] for x in orfs.values()]
        if aggregation_method == "lca-majority":
            res = tree.lca_star(taxonomies, majority_threshold=majority_threshold)
            contig_taxonomy = res["taxonomy"]
            error_function = res["pvalue"]
        elif aggregation_method == "lca":
            # TODO incorporate threshold into LCAs?
            contig_taxonomy = tree.lca(taxonomies)
            error_function = nettleton_pvalue(taxonomies, contig_taxonomy)
        # simple majority
        else:
            contig_taxonomy = BlastHits(taxonomies).majority()
            error_function = nettleton_pvalue(taxonomies, contig_taxonomy)
        lineage = ";".join(["%s" % (tree.tree[i].taxonomy) for i in tree.taxonomic_lineage(contig_taxonomy)])

        for idx in sorted(orfs.keys()):
            orf_function, orf_tax_id, bitscore, evalue = orfs[idx]
            orf_taxonomy = tree.tree[orf_tax_id].taxonomy
            print(contig, "%s_%s" % (contig, idx), lineage, error_function, orf_taxonomy,
                  orf_function, evalue, bitscore, sep="\t", file=output)


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


@cli.command("prepare-refseq", short_help="prepares refseq mapping files")
@click.argument("fasta", type=click.Path(exists=True))
@click.argument("namesdmp", type=click.Path(exists=True))
@click.argument("nodesdmp", type=click.Path(exists=True))
@click.argument("namemap", type=click.File("w"))
@click.argument("tree", type=click.File("w"))
def prepare_refseq_reference(fasta, namesdmp, nodesdmp, namemap, tree):
    """Takes the RefSeq FASTA which will look like:

        \b
        >gi|507053311|ref|WP_016124266.1| two-component sensor histidine kinase [Bacillus cereus]
        MNFNKRLIIQFIMQHVFVLVTLLIAVVAAFTYLIFLLTSTLYEPNIPDSDSFTISRYISSEDGHISLQSEVQDLIKEKND

    And from NCBI's taxonomy dump, names.dmp, which starts like:

        \b
        1 | all      |                         | synonym         |
        1 | root     |                         | scientific name |
        2 | Bacteria | Bacteria <prokaryotes>  | scientific name |

    And nodes.dmp:

        \b
        1 | 1      | no rank      |...
        2 | 131567 | superkingdom |...
        6 | 335928 | genus        |...

    To create a TSV map (NAMEMAP) of FASTA entry ID to function, taxonomy name, taxonomy ID, and parent
    taxonomy ID:

        \b
        gi|507053311|ref|WP_016124266.1|<tab>two-component sensor histidine kinase<tab>Bacillus cereus

    As well as TREE:

        \b
        Bacillus cereus<tab>1396<tab>86661

    """

    # most names, for translating reference fasta from tax name to tax id
    name_to_tax = {}
    tax_to_name = {}
    # for best names, for translating nodes.dmp tax id to tax name
    tax_to_scientific_name = {}

    logging.info("Reading in %s" % namesdmp)
    with open(namesdmp) as dmp:
        for tax_id, group in groupby(dmp, key=lambda x: [i.strip() for i in x.strip().split("|")][0]):
            scientific_name = ""
            synonym = ""
            for line in group:
                toks = [x.strip() for x in line.strip().split("|")]
                if toks[-2] == "scientific name":
                    tax_to_scientific_name[toks[0]] = toks[1]
                    break
                elif toks[-2] == "misspelling":
                    continue
                elif toks[-2] == "synonym":
                    synonym = toks[1]
                else:
                    if synonym:
                        tax_to_name[toks[0]] = synonym
                    else:
                        tax_to_name[toks[0]] = toks[1]
                name_to_tax[toks[1]] = toks[0]

    logging.info("Iterating over %s" % fasta)
    with open(fasta) as fa:
        for name, seq in read_fasta(fa):
            name_parts = name.partition(" ")
            # biotin--[acetyl-CoA-carboxylase] synthetase [Aeromonas allosaccharophila]
            function = name_parts[2].rpartition("[")[0].strip()
            taxonomy_name = name_parts[2].rpartition("[")[2].strip("]")
            # phosphoribosyltransferase [[Haemophilus] parasuis]
            if "[[" in name_parts[2]:
                taxonomy_name = "[%s" % taxonomy_name
            # https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?lvl=0&id=1737424
            if taxonomy_name == "Blautia sp. GD8":
                taxonomy_name = "Blautia sp. GD9"
            try:
                taxonomy_id = name_to_tax[taxonomy_name]
            except KeyError:
                taxonomy_name = taxonomy_name.replace("]", "").replace("[", "")
                taxonomy_id = name_to_tax[taxonomy_name]

            print(name_parts[0], function, taxonomy_id, sep="\t", file=namemap)

    logging.info("Reading in %s" % nodesdmp)
    with open(nodesdmp) as dmp:
        for line in dmp:
            toks = [x.strip() for x in line.strip().split("|")]
            try:
                taxonomy_name = tax_to_scientific_name[toks[0]]
            except KeyError:
                taxonomy_name = tax_to_name[toks[0]]
            print(toks[0], taxonomy_name, toks[1], sep="\t", file=tree)

    logging.info("Process complete")


@cli.command("prepare-eggnog", short_help="prepares eggnog mapping and reference files")
@click.argument("fasta", type=click.File("r"))
@click.argument("namemap", type=click.File("r"))
@click.argument("outputfasta", type=click.File("w"))
@click.argument("outputmap", type=click.File("w"))
def prepare_eggnog_reference(fasta, namemap, outputfasta, outputmap):
    names = set()
    # IDs for which we have uniprot AC
    for line in namemap:
        toks = line.strip().split("\t")
        names.add(toks[3])
    fasta_names = set()
    # IDs for which we have sequences
    for name, seq in read_fasta(fasta):
        name = name.partition(".")[-1]
        if name in names: continue
        fasta_names.add(name)
        print_fasta_record(name, seq, outputfasta)
    # the union of IDs
    for line in namemap:
        toks = line.strip().split("\t")
        if toks[3] in fasta_names:
            print(*toks, sep="\t", file=outputmap)


def print_fasta_record(name, seq, out_handle=sys.stdout, wrap=100):
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
    for header, group in groupby(fh, lambda line: line[0] == '>'):
        if header:
            line = next(group)
            name = line[1:].strip()
        else:
            seq = ''.join(line.strip() for line in group)
            yield name, seq


@cli.command("merge-tables", short_help="for a sample, merge its eggnog and refseq tables")
@click.argument("refseq", type=click.File("r"))
@click.argument("eggnog", type=click.File("r"))
@click.argument("output", type=click.File("w"))
def merge_tables(refseq, eggnog, output):
    """Takes the output from `refseq` and `eggnog` and combines them into a single TSV table.

    Headers are required and should contain 'contig' and 'orf' column labels.
    """
    import pandas as pd

    index_cols = ["contig", "orf"]
    try:
        ref_df = pd.read_table(refseq, index_col=index_cols)
    except ValueError:
        logging.critical("The expected headers ('contig', 'orf') are missing from %s" % refseq.name)
        sys.exit(1)
    logging.info("%d contained in %s" % (len(ref_df), refseq.name))
    try:
        egg_df = pd.read_table(eggnog, index_col=index_cols)
    except ValueError:
        logging.critical("The expected headers ('contig', 'orf') are missing from %s" % eggnog.name)
        sys.exit(1)
    logging.info("%d contained in %s" % (len(egg_df), eggnog.name))
    merged = pd.merge(left=ref_df, right=egg_df, how="outer", left_index=True, right_index=True)
    logging.info("%d total lines after merging" % len(merged))
    merged.to_csv(output, sep="\t", na_rep="NA")


if __name__ == "__main__":
    cli()
