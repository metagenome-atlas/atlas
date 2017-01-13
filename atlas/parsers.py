import contextlib
import gzip
import logging
import sqlite3
from atlas import BLAST6
from atlas.blast import BlastHits, Node, Tree, parse_blast_results_with_tree, process_orfs_with_tree, get_hit_from_blast_group
from atlas.utils import gzopen
from itertools import groupby


def cog_parser(tsv, namemap, output, summary_method='best', min_identity=60, min_bitscore=0, min_length=60, max_evalue=0.000001, top_fraction=1, max_hits=10, table_name="cog"):
    """Parse BLAST hits from COG reference database.

    The BLAST hits are assumed to be sorted by query with decreasing bitscores (best alignment first):

        \b
        sort -k1,1 -k12,12rn tsv > sorted_tsv

    Expected columns in the COG database:

        \b
        cog_protein_key (unique and matches sequence name in reference fasta)
        cog_protein_id
        cog_id
        cog_functional_class
        cog_annotation
        cog_functional_class_description

    Args:
        tsv (str): blast hits file path in default tabular format
        namemap (str): sqlite database file path
        output (str): :py:class:`click.File` in write mode
        summary_method (str): either 'majority' or 'best'; summary method for annotating ORFs; when majority and there is no majority, best is used
        min_identity (int): minimum allowable percent ID of BLAST hit
        min_bitscore (int): minimum allowable bitscore of BLAST hit; 0 disables
        min_length (int): minimum allowable BLAST alignment length
        max_evalue (float): maximum allowable e-value of BLAST hit
        top_fraction (float): filters ORF BLAST hits before finding majority by only keep hits within this fraction, e.g. 0.98, of the highest bitscore
        max_hits (int): maximum number of BLAST hits to consider when summarizing ORFs as a majority
        table_name (str): table name within namemap database; expected columns are listed above

    """
    logging.info("Parsing %s" % tsv)

    if top_fraction == 1:
        top_fraction = None

    print("contig", "orf", "cog_protein_id", "cog_id", "cog_functional_class", "cog_annotation",
          "cog_functional_class_description", "%s_evalue" % table_name, "%s_bitscore" % table_name,
          sep="\t", file=output)
    with contextlib.closing(sqlite3.connect(namemap)) as conn, gzopen(tsv) as blast_tab_fh:
        cursor = conn.cursor()
        for query, qgroup in groupby(blast_tab_fh, key=lambda x: x.partition("\t")[0]):

            contig_name, _, orf_idx = query.rpartition("_")
            hit_id, evalue, bitscore = get_hit_from_blast_group(qgroup, max_hits, top_fraction,
                                                                min_length, min_identity,
                                                                max_evalue, min_bitscore,
                                                                summary_method)
            cog_protein_id = "NA"
            cog_id = "NA"
            cog_functional_class = "NA"
            cog_annotation = "NA"
            cog_functional_class_description = "NA"

            # everything could have been filtered out due to user constraints
            if hit_id:
                cursor.execute('SELECT cog_protein_id, cog_id, cog_functional_class, \
                                       cog_annotation, cog_functional_class_description \
                                FROM %s \
                                WHERE cog_protein_key="%s"' % (table_name, hit_id))
                cog_protein_id, cog_id, cog_functional_class, cog_annotation, \
                    cog_functional_class_description = cursor.fetchone()
            print(contig_name, "%s_%s" % (contig_name, orf_idx), cog_protein_id, cog_id,
                  cog_functional_class, cog_annotation, cog_functional_class_description,
                  evalue, bitscore, sep="\t", file=output)
    logging.info("Complete")


def cazy_parser(tsv, namemap, output, summary_method='best', min_identity=60, min_bitscore=0, min_length=60, max_evalue=0.000001, top_fraction=1, max_hits=10, table_name="cazy"):
    """Parse BLAST hits from CAZy reference database.

    The BLAST hits are assumed to be sorted by query with decreasing bitscores (best alignment first):

        \b
        sort -k1,1 -k12,12rn tsv > sorted_tsv

    Expected columns in the CAZy database:

        \b
        cazy_gene
        cazy_family
        cazy_class
        cazy_ec

    For a given gene match, all possible ECs are returned in a single line separated by '|'.

    Args:
        tsv (str): blast hits file path in default tabular format
        namemap (str): sqlite database file path
        output (str): :py:class:`click.File` in write mode
        summary_method (str): either 'majority' or 'best'; summary method for annotating ORFs; when majority and there is no majority, best is used
        min_identity (int): minimum allowable percent ID of BLAST hit
        min_bitscore (int): minimum allowable bitscore of BLAST hit; 0 disables
        min_length (int): minimum allowable BLAST alignment length
        max_evalue (float): maximum allowable e-value of BLAST hit
        top_fraction (float): filters ORF BLAST hits before finding majority by only keep hits within this fraction, e.g. 0.98, of the highest bitscore
        max_hits (int): maximum number of BLAST hits to consider when summarizing ORFs as a majority
        table_name (str): table name within namemap database; expected columns are listed above

    """
    logging.info("Parsing %s" % tsv)

    if top_fraction == 1:
        top_fraction = None

    print("contig", "orf", "cazy_gene", "cazy_family", "cazy_class", "cazy_ec",
          "%s_evalue" % table_name, "%s_bitscore" % table_name, sep="\t", file=output)
    with contextlib.closing(sqlite3.connect(namemap)) as conn, gzopen(tsv) as blast_tab_fh:
        cursor = conn.cursor()
        for query, qgroup in groupby(blast_tab_fh, key=lambda x: x.partition("\t")[0]):

            contig_name, _, orf_idx = query.rpartition("_")
            hit_id, evalue, bitscore = get_hit_from_blast_group(qgroup, max_hits, top_fraction, min_length, min_identity, max_evalue, min_bitscore, summary_method)
            cazy_gene = "NA"
            cazy_family = "NA"
            cazy_class = "NA"
            cazy_ec = "NA"

            # everything could have been filtered out due to user constraints
            if hit_id:
                cursor.execute('SELECT cazy_gene, cazy_family, cazy_class, cazy_ec \
                                FROM %s \
                                WHERE cazy_gene="%s"' % (table_name, hit_id))
                cazy_gene, cazy_family, cazy_class, cazy_ec = cursor.fetchone()
            print(contig_name, "%s_%s" % (contig_name, orf_idx), cazy_gene, cazy_family,
                  cazy_class, cazy_ec, evalue, bitscore, sep="\t", file=output)
    logging.info("Complete")


def eggnog_parser(tsv, namemap, output, summary_method, min_identity, min_bitscore, min_length, max_evalue, top_fraction, max_hits, table_name):
    """Parse BLAST hits from EGGNOG.

    The BLAST hits are assumed to be sorted by query with decreasing bitscores (best alignment first):

        \b
        sort -k1,1 -k12,12rn tsv > sorted_tsv

    Expected columns in the EggNOG database:

        \b
        uniprot_ac
        eggnog_ssid_b
        eggnog_species_id
        uniprot_id
        ko_id
        ko_level1_name
        ko_level2_name
        ko_level3_id
        ko_level3_name
        ko_gene_symbol
        ko_product
        ko_ec

    Args:
        tsv (str): blast hits file path in default tabular format
        namemap (str): sqlite database file path
        output (str): :py:class:`click.File` in write mode
        summary_method (str): either 'majority' or 'best'; summary method for annotating ORFs; when majority and there is no majority, best is used
        min_identity (int): minimum allowable percent ID of BLAST hit
        min_bitscore (int): minimum allowable bitscore of BLAST hit; 0 disables
        min_length (int): minimum allowable BLAST alignment length
        max_evalue (float): maximum allowable e-value of BLAST hit
        top_fraction (float): filters ORF BLAST hits before finding majority by only keep hits within this fraction, e.g. 0.98, of the highest bitscore
        max_hits (int): maximum number of BLAST hits to consider when summarizing ORFs as a majority
        table_name (str): table name within namemap database; expected columns are listed above

    """
    logging.info("Parsing %s" % tsv)

    if top_fraction == 1:
        top_fraction = None

    print("contig", "orf", "uniprot_ac", "eggnog_ssid_b", "eggnog_species_id", "uniprot_id",
          "ko_id", "ko_level1_name", "ko_level2_name", "ko_level3_id", "ko_level3_name",
          "ko_gene_symbol", "ko_product", "ko_ec", "%s_evalue" % table_name,
          "%s_bitscore" % table_name, sep="\t", file=output)

    with contextlib.closing(sqlite3.connect(namemap)) as conn, gzopen(tsv) as blast_tab_fh:
        cursor = conn.cursor()
        for query, qgroup in groupby(blast_tab_fh, key=lambda x: x.partition("\t")[0]):

            contig_name, _, orf_idx = query.rpartition("_")
            hit_id, evalue, bitscore = get_hit_from_blast_group(qgroup, max_hits, top_fraction, min_length, min_identity, max_evalue, min_bitscore, summary_method)
            uniprot_ac = "NA"
            eggnog_ssid_b = "NA"
            eggnog_species_id = "NA"
            uniprot_id = "NA"
            ko_id = "NA"
            ko_level1_name = "NA"
            ko_level2_name = "NA"
            ko_level3_id = "NA"
            ko_level3_name = "NA"
            ko_gene_symbol = "NA"
            ko_product = "NA"
            ko_ec = "NA"

            # everything could have been filtered out due to user constraints
            if hit_id:
                cursor.execute('SELECT uniprot_ac, eggnog_ssid_b, eggnog_species_id, uniprot_id, \
                                       ko_id, ko_level1_name, ko_level2_name, \
                                       ko_level3_id, ko_level3_name, ko_gene_symbol, ko_product, ko_ec \
                                FROM %s \
                                WHERE eggnog_ssid_b="%s"' % (table_name, hit_id))
                try:
                    uniprot_ac, eggnog_ssid_b, eggnog_species_id, uniprot_id, ko_id, \
                        ko_level1_name, ko_level2_name, ko_level3_id, ko_level3_name, \
                        ko_gene_symbol, ko_product, ko_ec = cursor.fetchone()
                # legacy before database was pruned; can have hits not in metadata
                except TypeError:
                    logging.warning("'%s' not present in database" % hit_id)
                    pass

            # print for this query
            print(contig_name, "%s_%s" % (contig_name, orf_idx), uniprot_ac, eggnog_ssid_b,
                  eggnog_species_id, uniprot_id, ko_id, ko_level1_name, ko_level2_name,
                  ko_level3_id, ko_level3_name, ko_gene_symbol, ko_product, ko_ec, evalue,
                  bitscore, sep="\t", file=output)
    logging.info("Complete")


def expazy_parser(tsv, namemap, output, summary_method, min_identity, min_bitscore, min_length,
                   max_evalue, top_fraction, max_hits, table_name):
    """Parse BLAST hits from ExPAZy reference database.

    The BLAST hits are assumed to be sorted by query with decreasing bitscores (best alignment first):

        \b
        sort -k1,1 -k12,12rn tsv > sorted_tsv

    Expected columns in the ExPAZy database:

        \b
        uniprot_entry
        uniparc_entry
        expazy_ec
        expazy_name

    For a given UniParc match, all possible ECs and ExPAZy recommended names are returned in a
    single line separated by '|'.

    Args:
        tsv (str): blast hits file path in default tabular format
        namemap (str): sqlite database file path
        output (str): :py:class:`click.File` in write mode
        summary_method (str): either 'majority' or 'best'; summary method for annotating ORFs; when majority and there is no majority, best is used
        min_identity (int): minimum allowable percent ID of BLAST hit
        min_bitscore (int): minimum allowable bitscore of BLAST hit; 0 disables
        min_length (int): minimum allowable BLAST alignment length
        max_evalue (float): maximum allowable e-value of BLAST hit
        top_fraction (float): filters ORF BLAST hits before finding majority by only keep hits within this fraction, e.g. 0.98, of the highest bitscore
        max_hits (int): maximum number of BLAST hits to consider when summarizing ORFs as a majority
        table_name (str): table name within namemap database; expected columns are listed above

    """
    logging.info("Parsing %s" % tsv)

    if top_fraction == 1:
        top_fraction = None

    print("contig", "orf", "expazy_ec", "expazy_name", "%s_evalue" % table_name,
          "%s_bitscore" % table_name, sep="\t", file=output)
    with contextlib.closing(sqlite3.connect(namemap)) as conn, gzopen(tsv) as blast_tab_fh:
        cursor = conn.cursor()
        for query, qgroup in groupby(blast_tab_fh, key=lambda x: x.partition("\t")[0]):

            contig_name, _, orf_idx = query.rpartition("_")
            hit_id, evalue, bitscore = get_hit_from_blast_group(qgroup, max_hits, top_fraction, min_length, min_identity, max_evalue, min_bitscore, summary_method)
            ecs = "NA"
            names = "NA"

            # everything could have been filtered out due to user constraints
            if hit_id:
                cursor.execute('SELECT expazy_ec, expazy_name \
                                FROM %s \
                                WHERE uniparc_entry="%s"' % (table_name, hit_id))
                ecs, names = cursor.fetchone()
            print(contig_name, "%s_%s" % (contig_name, orf_idx), ecs, names, evalue, bitscore, sep="\t", file=output)
    logging.info("Complete")


def refseq_parser(tsv, namemap, treefile, output, summary_method, aggregation_method, majority_threshold, min_identity, min_bitscore, min_length, max_evalue, max_hits, table_name, top_fraction):
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
    orf_assignments = parse_blast_results_with_tree(tsv, namemap, summary_method=summary_method, tree=tree,
                          min_identity=min_identity, min_bitscore=min_bitscore,
                          min_length=min_length, max_evalue=max_evalue,
                          max_hits_per_orf=max_hits, table_name=table_name,
                          top_fraction_of_hits=top_fraction)
    logging.info("Assigning taxonomies to contigs using %s" % aggregation_method)
    process_orfs_with_tree(orf_assignments, tree, output, aggregation_method, majority_threshold, table_name)
    logging.info("Complete")
