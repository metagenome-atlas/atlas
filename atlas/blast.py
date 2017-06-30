import bisect
import contextlib
import logging
import sqlite3
from atlas import BLAST6, TAX_LEVELS
from atlas.utils import gzopen, index_of_list_items, nettleton_pvalue
from collections import Counter, defaultdict, deque, OrderedDict
from itertools import groupby


class Node(object):

    def __init__(self, taxonomy, node_id, parent_id, tax_level):
        """Represents a node within a tree.

        Args:
            taxonomy (str): taxonomy name or ID
            node_id (str): taxonomy ID
            parent_id (str): taxonomy ID of parent
            tax_level (str): the taxonomic level for this node_id

        """
        # the current node's string ID
        self.taxonomy = taxonomy
        # the current node's digit ID
        self.node_id = node_id
        self.parent_id = parent_id
        self.tax_level = tax_level


class Tree(object):

    def __init__(self, tree_file):
        """Builds reference dictionary of Taxonomy Name, Taxonomy ID, and Parent Taxonomy ID."""
        self.tree = defaultdict(dict)

        with open(tree_file) as tree_fh:
            for line in tree_fh:
                toks = line.strip().split("\t")
                if not toks[0] == '1' and not toks[2] == '1':
                    assert not toks[0] == toks[2]
                if not len(toks) == 4:
                    logging.warning("Line [%s] does not have ID, NAME, PARENTID, TAX LEVEL" % line.strip())
                    continue
                self.add_node(toks[1], toks[0], toks[2], toks[3])

    def add_node(self, taxonomy, node_id, parent_id, tax_level):
        """Adds node to tree dictionary.

        Args:
            taxonomy (str): the taxonomy name
            node_id (str): the taxonomy id
            parent_id (str): the parent's taxonomy id
            tax_level (str): the taxonomic level for this node_id

        """
        # taxonomy id to node mapping; ensures unique nodes despite non-unique names
        self.tree[node_id] = Node(taxonomy, node_id, parent_id, tax_level)

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


def parse_blast_results_with_tree(blast_tab, name_map, summary_method, tree, min_identity=70,
                                  min_bitscore=0, min_length=60, max_evalue=0.000001,
                                  max_hits_per_orf=10, top_fraction_of_hits=None,
                                  table_name="refseq", lca_threshold=1):
    """Parse BLAST results (-outfmt 6), filter, and aggregate ORF taxonomies.

    Args:
        blast_tab (str): file path to blast TSV file
        name_map (dict): dict of tuples from parse_tree_annotation
        summary_method (dict): method of ORF annotation selection

        lca_threshold (float): the first parent above this fraction of representation (its count is
            greater than the total * lca_threshold)

    Returns:
        dict: dict of dicts where first key is contig name, inner key is ORF ID; values are tuple of
              protein function, taxonomy ID, bitscore, evalue

    Raises:
        AssertionError when ORF summary method is not supported (['lca', 'best', 'majority'])
    """

    # allowing 1 and 0 to disable
    if top_fraction_of_hits == 1:
        top_fraction_of_hits = None

    assert summary_method in ["lca", "best", "majority"]
    contigs = defaultdict(dict)

    with contextlib.closing(sqlite3.connect(name_map)) as conn, gzopen(blast_tab) as blast_tab_fh:
        cursor = conn.cursor()
        # group hits by ORF (column 2)
        for orf_id, qgroup in groupby(blast_tab_fh, key=lambda x: x.split("\t")[1]):

            protein_function = "hypothetical protein"
            protein_set = False
            taxonomy_id = "1"
            bitscore = "NA"
            evalue = "NA"
            orf_hits = BlastHits(max_hits=max_hits_per_orf, top_fraction=top_fraction_of_hits)
            lines = []

            # iterate over blast hits per ORF
            for hsp in qgroup:
                # HSPs will now have contig in column 1
                toks = hsp.strip().split("\t")
                # remove extra column from toks
                contig_name = toks.pop(0)
                # convert toks to dictionary
                toks = dict(zip(BLAST6, toks))

                if (int(toks["length"]) < min_length or
                        float(toks["pident"]) < min_identity or
                        float(toks["evalue"]) > max_evalue):
                    continue
                if min_bitscore and float(toks["bitscore"]) < min_bitscore:
                    # input is sorted by decreasing bitscore
                    break

                cursor.execute('SELECT function, taxonomy FROM %s WHERE name="%s"' % (table_name, toks["sseqid"]))
                current_function, current_taxonomy = cursor.fetchone()

                # update taxonomy based on pident; would be similar to 16S taxonomy assignments
                # current_taxonomy = tree.climb_tree(current_taxonomy, float(toks["pident"]))

                if summary_method == "best":
                    taxonomy_id = current_taxonomy
                    protein_function = current_function
                    bitscore = toks["bitscore"]
                    evalue = toks["evalue"]
                    break

                # TODO implement bitscore ratio as a measure of alignment quality as a function of input sequence
                orf_hits.add(current_taxonomy, toks["bitscore"])
                toks["current_function"] = current_function
                toks["current_taxonomy"] = current_taxonomy
                lines.append(toks)

            # summary method is majority and we have passing HSPs
            if not summary_method == "best" and lines:
                if summary_method == "majority":
                    taxonomy_id = orf_hits.majority()
                    for toks in lines:
                        if toks["current_taxonomy"] == taxonomy_id:
                            bitscore = toks["bitscore"]
                            evalue = toks["evalue"]
                            protein_function = toks["current_function"]
                            break
                # summary method is 'lca'
                else:
                    orf_hits.names.reverse()
                    taxonomy_id = tree.lca(orf_hits.names, threshold=lca_threshold)
                    # grabbing best hit's bitscore and evalue
                    bitscore = lines[0]["bitscore"]
                    evalue = lines[0]["evalue"]
                    protein_function = lines[0]["current_function"]

                if bitscore == "NA":
                    logging.critical("The summarized ID (%s) was not assigned metadata" % taxonomy_id)

            contigs[contig_name][orf_id] = (protein_function, taxonomy_id, bitscore, evalue)

    return contigs


def validate_lineage(lineage, sep=";"):
    """
    >>> lineage = {"p":"Basidiomycota","c":"Tremellomycetes","o":"Tremellales","g":"Cryptococcus"}
    >>> validate_lineage(lineage)
    'k__?;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__?;g__Cryptococcus;s__?'
    """
    levels = ["k" if tax_level == "superkingdom" else tax_level[0] for tax_level in TAX_LEVELS]
    valid_lineage = []
    for idx in levels:
        # removes commas in tax names
        valid_lineage.append("%s__%s" % (idx, lineage.get(idx, "?").replace(",", "")))
    return sep.join(valid_lineage)


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
        lineage = {}
        for item in tree.taxonomic_lineage(contig_taxonomy):
            node = tree.tree[item]
            if node.tax_level in TAX_LEVELS:
                # does not account for "no rank" and some other cases of "unclassified"
                lineage["k" if node.tax_level == "superkingdom" else node.tax_level[0]] = node.taxonomy
        lineage = validate_lineage(lineage)

        for idx in sorted(orfs.keys()):
            orf_function, orf_tax_id, bitscore, evalue = orfs[idx]
            orf_taxonomy = tree.tree[orf_tax_id].taxonomy
            print(contig, idx, lineage, error_function, orf_taxonomy,
                  orf_function, evalue, bitscore, sep="\t", file=output)
