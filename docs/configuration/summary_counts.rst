Summary Counts
==============

ATLAS will generate a single annotation table with all annotation metadata
and counts per open reading frame. ``summary_counts`` is configurable to
map count data across the annotation metadata in various ways.

To summarize ENZYME ECs across species, you will have needed to use RefSeq and
ENZYME databases then use::

    summary_counts:
        taxonomy:
            levels:
                - species
            ENZYME:
                - enzyme_name
                - enzyme_ec

Taxonomic levels include:

    + kingdom
    + domain
    + phylum
    + class
    + order
    + family
    + genus
    + species

Another example is when you do not need taxonomy and just want to see a
breakdown across CAZy classes::

    summary_counts:
        CAZy_Families:
            - cazy_family
            - cazy_class

Multiple tables can be specified for taxonomic levels and for non-taxonomic
annotations::

    summary_counts:
        taxonomy:
            levels:
                - phylum
                - class
                - order
                - species
            COG:
                - cog_id
            CAZy_family:
                - cazy_family
            ENZYME:
                - enzyme_name
                - enzyme_ec
        CAZY_EC:
            - cazy_ec
        COG:
            - cog_id
            - cog_product
        ENZYME:
            - enzyme_name
            - enzyme_ec
