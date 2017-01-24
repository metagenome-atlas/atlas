Configuration
=============




## Annotation Parameters

annotation:
    ## ORFs
    # https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    translation_table: 11
    # when counting reads aligning to ORFs, require at least this many bp
    # overlapping the ORF
    minimum_overlap: 20

    references:
        eggnog:
            # non-tree based reference requires namemap database and fasta
            namemap: /pic/projects/mint/atlas_databases/functional/eggnog/eggnog4_nonredundant.db
            fasta: /pic/projects/mint/atlas_databases/functional/eggnog/eggnog4_nonredundant.fasta
            # number of entries per FASTA to be aligned with DIAMOND
            chunk_size: 250000
            # 'fast' or 'sensitive'
            run_mode: fast
            # setting top_seqs to 5 will report all alignments whose score is
            # at most 5% lower than the top alignment score for a query
            top_seqs: 5
            # maximum e-value to report alignments
            e_value: "0.000001"
            # minimum identity % to report an alignment
            min_identity: 50
            # minimum query cover % to report an alignment
            query_coverage: 60
            # gap open penalty
            gap_open: 11
            # gap extension penalty
            gap_extend: 1
            # Block size in billions of sequence letters to be processed at a time.
            # This is the main parameter for controlling DIAMOND's memory usage.
            # Bigger numbers will increase the use of memory and temporary disk space,
            # but also improve performance. The program can be expected to roughly use
            # six times this number of memory (in GB).
            block_size: 4
            # The number of chunks for processing the seed index (default=4). This
            # option can be additionally used to tune the performance. It is
            # recommended to set this to 1 on a high memory server, which will
            # increase performance and memory usage, but not the usage of temporary
            # disk space.
            index_chunks: 4
            # 'majority' or 'best'; summary method for annotating ORFs
            summary_method: best
            # minimum allowable BLAST alignment length
            min_length: 60
            # maximum allowable e-value of BLAST hit when parsing DIAMOND hits
            max_evalue: 0.000001
            # maximum number of BLAST hits to consider when summarizing ORFs
            max_hits: 10
            # filters ORF BLAST hits by only keep hits within this fraction of
            # the highest bitscore; this is recommended over max_hits
            top_fraction: 0.50
            # minimum allowable BLAST alignment bitscore; 0 effectively disables
            min_bitscore: 0
        refseq:
            # tree based reference requires namemap database, tree, and fasta
            namemap: /pic/projects/mint/atlas_databases/functional/refseq/refseq78.complete.nonredundant_protein.faa.db
            tree: /pic/projects/mint/atlas_databases/functional/refseq/refseq78.complete.nonredundant_protein.faa.tree
            fasta: /pic/projects/mint/atlas_databases/functional/refseq/refseq78.complete.nonredundant_protein.faa
            # number of entries per FASTA to be aligned with DIAMOND
            chunk_size: 250000
            run_mode: fast
            top_seqs: 5
            e_value: "0.000001"
            min_identity: 50
            query_coverage: 60
            gap_open: 11
            gap_extend: 1
            block_size: 6
            index_chunks: 1
            # 'lca', 'majority', or 'best'; summary method for annotating ORFs; when
            # using LCA, it's recommended that one limits the number of hits using a
            # low top_fraction
            summary_method: best
            # 'lca', 'lca-majority', or 'majority'; summary method for aggregating ORF
            # taxonomic assignments to contig level assignment; 'lca' will result in
            # most stringent, least specific assignments
            aggregation_method: lca-majority
            # constitutes a majority fraction at tree node for 'lca-majority' ORF
            # aggregation method
            majority_threshold: 0.51
            # minimum allowable BLAST alignment length
            min_length: 60
            # maximum allowable e-value of BLAST hit
            max_evalue: 0.000001
            # maximum number of BLAST hits to consider when summarizing ORFs; can
            # drastically alter ORF LCA assignments if too high without further limits
            max_hits: 10
            top_fraction: 0.50
        expazy:
            namemap: /pic/projects/mint/atlas_databases/functional/expazy/expazy.db
            fasta: /pic/projects/mint/atlas_databases/functional/expazy/expazy.fasta
            chunk_size: 500000
            # 'fast' or 'sensitive'
            run_mode: fast
            top_seqs: 2
            index_chunks: 1
            summary_method: majority
        cazy:
            namemap: /pic/projects/mint/atlas_databases/functional/dbcan/dbcan.db
            fasta: /pic/projects/mint/atlas_databases/functional/dbcan/dbcan.fasta
            chunk_size: 500000
            # 'fast' or 'sensitive'
            run_mode: fast
            top_seqs: 2
            index_chunks: 1
            summary_method: majority

summary_counts:
    # Possible columns table column values upon which to aggregate:
        # contig, orf

        # from refseq:
        # taxonomy, orf_taxonomy, refseq_product

        # from eggnog:
        # uniprot_ac, eggnog_ssid_b, eggnog_species_id, uniprot_id, cog_func_id, cog_id,
        # cog_product, cog_level1_code, cog_level1_name, cog_level2_name,
        # ko_id, ko_level1_name, ko_level2_name, ko_level3_id,
        # ko_level3_name, ko_gene_symbol, ko_product, ko_ec

        # from expazy:
        # expazy_name, expazy_ec

        # from cazy (dbcan):
        # cazy_gene, cazy_family, cazy_class, cazy_ec

    # this is a special case to allow for taxon level specification
    taxonomy:
        # limit taxonomy in classification to the depth specified
        # possible values: kingdom, domain, phylum, class, order, family, genus, species
        # all levels if omitted
        levels:
            - phylum
            - class
            - order
            - species
        # tables to generate at these taxonomic levels
        KO:
            - ko_id
            - ko_ec
        COG:
            - cog_id
        CAZy_EC:
            - cazy_ec
        CAZy_family:
            - cazy_family
        ExPAZy:
            - expazy_name
            - expazy_ec
    KO:
        - ko_id
        - ko_gene_symbol
        - ko_product
        - ko_ec
    KO_lvl1:
        - ko_level1_name
    KO_lvl2:
        - ko_level2_name
    KO_lvl3:
        - ko_level3_name
    CAZY_EC:
        - cazy_ec
    COG:
        - cog_id
        - cog_product
    COG_lvl1:
        - cog_level1_name
        - cog_level2_name
    ExPAZy:
        - expazy_name
        - expazy_ec


# Output
