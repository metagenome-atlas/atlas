Summary Counts
==============

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

        # from ENZYME:
        # enzyme_name, enzyme_ec

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
        ENZYME:
            - enzyme_name
            - enzyme_ec
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
    ENZYME:
        - enzyme_name
        - enzyme_ec


# Output
