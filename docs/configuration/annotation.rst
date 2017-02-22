Annotation
==========

Translation Table
-----------------

By default, translation table 11 is used to find open reading frames among
passing contig sequences. Other codes are available at
https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi.


Minimum Overlap
---------------

When counting reads overlapping coding sequence, require this much read
overlap.

**Default: 1**

::

    annotation:
        minimum_overlap: 1


Restricting Counts
------------------

Counts can be restricted to primary alignments only using ``primary_only`` in
addition to being able to control behavior associated with multi-mapped reads.
As the alignment stage does allow up to 10 possible hits per sequence, one may
want to restrict their counts using ``multi_mapping``.

Defaults are reflected in the example::

    annotation:
        primary_only: false
        multi_mapping: true


``MaxBin`` Options
------------------

Used to set iterations, minimum contig length to be considered, and the
probability threshold for EM final classification.

::

    annotation:
        maxbin_max_iteration: 50
        maxbin_min_contig_length: 500
        maxbin_prob_threshold: 0.9



References
----------

Any combination of these annotation tables works, but each provides data that
can later be subset downstream. Any column later specified in summary tables
requires that its reference be specified. For instance, specifying taxonomic
breakdowns across EC will require RefSeq for the taxonomic assignments and
ENZYME for the most specific EC data.

For each reference, one can specify options that are sent to DIAMOND and others
that are used when selecting the best annotation hit or performing an LCA.


EggNOG version 4.5
``````````````````

EggNOG output columns:

+----------------+--------------------------+
| Value          | Definition               |
+================+==========================+
| uniprot_ac     | UniProt accession number |
+----------------+--------------------------+
| eggnog_ssid_b  | EggNOG accession         |
+----------------+--------------------------+
| uniprot_id     | UniProt name             |
+----------------+--------------------------+
| ko_id          | KO                       |
+----------------+--------------------------+
| ko_level1_name | KEGG level 1 name        |
+----------------+--------------------------+
| ko_level2_name | KEGG level 2 name        |
+----------------+--------------------------+
| ko_level3_id   | KEGG pathway ID          |
+----------------+--------------------------+
| ko_level3_name | KEGG pathway name        |
+----------------+--------------------------+
| ko_gene_symbol | Gene symbol              |
+----------------+--------------------------+
| ko_product     | Protein product          |
+----------------+--------------------------+
| ko_ec          | EC                       |
+----------------+--------------------------+

::

    annotation:
        references:
            eggnog:
                namemap: /path/to/eggnog.db
                dmnd: /path/to/eggnog.dmnd


Nonredundant RefSeq version 78
``````````````````````````````

RefSeq is the largest reference provided in ATLAS and is used to assign product
and taxonomy to open reading frames. Those ORF assignments are then used to
calculate a taxonomy of the contig.

The reference consists of a map, tree, and DIAMOND formatted database.

+----------------+------------------------------------+
| Value          | Definition                         |
+================+====================================+
| taxonomy       | RefSeq taxonomy assigned to contig |
+----------------+------------------------------------+
| orf_taxonomy   | RefSeq taxonomy assigned to ORF    |
+----------------+------------------------------------+
| refseq_product | RefSeq product of ORF              |
+----------------+------------------------------------+

The RefSeq database also requires an additional reference argument, `tree`,
in addition to `namemap` and `dmnd` when specifying database files.

::

    annotation:
        references:
            refseq:
                namemap: /path/to/refseq.db
                tree: /path/to/refseq.tree
                dmnd: /path/to/refseq.dmnd


ENZYME
``````

Last updated on January 30, 2017.

Part of the `ExPASy portal`_, ENZYME provides a
complete EC protein sequence reference.

.. _ExPASy portal: http://enzyme.expasy.org/

+-------------+-------------------+
| Value       | Definition        |
+=============+===================+
| enzyme_name | The accepted name |
+-------------+-------------------+
| enzyme_ec   | EC                |
+-------------+-------------------+

Enzymes will be assigned as complete ECs and not broken down by level, such
that no assignment will contain a dash.

::

    annotation:
        references:
            enzyme:
                namemap: /path/to/enzyme.db
                dmnd: /path/to/enzyme.dmnd


dbCAN version 5.0 (CAZy)
````````````````````````

Provides CAZy reference annotation via dbCAN_.

.. _dbCAN: http://csbl.bmb.uga.edu/dbCAN/

+-------------+-----------------------------------------+
| Value       | Definition                              |
+=============+=========================================+
| cazy_gene   | Associated gene name of origin sequence |
+-------------+-----------------------------------------+
| cazy_family | CAZy family name                        |
+-------------+-----------------------------------------+
| cazy_class  | CAZy class name                         |
+-------------+-----------------------------------------+
| cazy_ec     | Associated EC; may contain dashes       |
+-------------+-----------------------------------------+

::

    annotation:
        references:
            cazy:
                namemap: /path/to/cazy.db
                dmnd: /path/to/cazy.dmnd


COG
```

Version is COG2014 via COG_.

.. _COG: ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data

+----------------------------------+-------------------------------+
| Value                            | Definition                    |
+==================================+===============================+
| cog_protein_id                   | Protein ID of origin sequence |
+----------------------------------+-------------------------------+
| cog_id                           | COG ID, e.g. COG0593          |
+----------------------------------+-------------------------------+
| cog_functional_class             | Functional class code         |
+----------------------------------+-------------------------------+
| cog_annotation                   | Functional class name         |
+----------------------------------+-------------------------------+
| cog_functional_class_description | Functional class description  |
+----------------------------------+-------------------------------+

::

    annotation:
        references:
            cog:
                namemap: /path/to/cog.db
                dmnd: /path/to/cog.dmnd


Reference Options
-----------------

Within each reference database, the user has the flexibility to optimize
performance across their compute environment and control the number of
alignment hits in various ways.


Name Map
````````

A SQL database that maps the fasta reference sequence name to the reference
metadata.

::

    annotation:
        references:
            refseq:
                namemap: /path/to/refseq.db


Reference Taxonomy Tree
```````````````````````

This applies only to the RefSeq database at this point as it is the only
reference utilizing a tree.

::

    annotation:
        references:
            refseq:
                tree: /path/to/refseq.tree


Alignment Index
```````````````

A DIAMOND formatted alignment index per reference database.

::

    annotation:
        references:
            refseq:
                dmnd: /path/to/refseq.dmnd


Query Chunk Size
````````````````

The number of entries per FASTA to be aligned with DIAMOND. To maximize large
compute clusters, this can be reduced to 100k to split input queries across
multiple blades. 250k rarely splits contigs into multiple input files, but
if you desire no splitting, set to something higher like 500000.

**Default: 250000**

::

    annotation:
        references:
            refseq:
                chunk_size: 250000


Run Mode
````````

DIAMOND alignment mode. Either 'fast' of 'more-sensitive'.

**Default: fast**

::

    annotation:
        references:
            refseq:
                run_mode: fast


Top Percent of Sequences
````````````````````````

Applies to reported local alignments. One can later further filter hits by
setting summary specific metrics. This allows more downstream customization
options and re-running of the protocol without having to redo the
computationally expensive annotation step.

**Default: 5**

::

    annotation:
        references:
            refseq:
                top_seqs: 10


Alignment e-value
`````````````````

Similar to the above, this filter applies to only the alignment step and
filters based on e-value.

**Default: 0.000001**

::

    annotation:
        references:
            refseq:
                e_value: 0.000001


Minimum Identity
````````````````

Filters DIAMOND hits based on minimum matching identity percentage.

**Default: 50**

::

    annotation:
        references:
            refseq:
                min_identity: 50


Query Coverage
``````````````

Require this much of the query sequence to be matched above ``min_identity``.

**Default: 60**

::

    annotation:
        references:
            refseq:
                query_coverage: 60


Gap Open Penalty
````````````````

A lower gap open penalty may allow more possible, lower identity hits.

**Default: 11**

::

    annotation:
        references:
            refseq:
                gap_open: 11


Gap Extend Penalty
``````````````````

A higher extend penalty will reduce allowable indel lengths in matches.

**Default: 1**

::

    annotation:
        references:
            refseq:
                gap_extend: 1


Block Size
``````````

Block size in billions of sequence letters to be processed at a time.
This is the main parameter for controlling DIAMOND's memory usage.
Bigger numbers will increase the use of memory and temporary disk space,
but also improve performance. The program can be expected to roughly use
six times this number of memory (in GB).

**Default: 2**

::

    annotation:
        references:
            refseq:
                block_size: 1


Index Chunks
````````````

The number of chunks for processing the seed index. This option can be
additionally used to tune the performance. It is recommended to set this to 1
on a high memory server, which will increase performance and memory usage, but
not the usage of temporary disk space.

**Default: 4**

::

    annotation:
        references:
            refseq:
                index_chunks: 4


Summary Method
``````````````

This is the summary method for annotating open reading frames. 'lca' performs
an LCA on the hits which can be limited using ``max_hits`` and
``top_fraction``. Other options are 'majority' which takes the majority target
hit after filtering alignments and 'best' which simply chooses the top hit.

**Default: best**

::

    annotation:
        references:
            refseq:
                summary_method: best


Aggregation Method
``````````````````

The summary method for aggregating ORF taxonomic assignments to a contig level
assignment.

+--------------+--------------------------------------------------------------+
| Value        | Definition                                                   |
+==============+==============================================================+
| lca-majority | Taxonomy is based on counts at tree nodes and works in       |
|              | combination with ``majority_threshold``; ``lca-majority`` is |
|              | a balance between the very stringent ``lca`` and the least   |
|              | restrictive ``majority``                                     |
+--------------+--------------------------------------------------------------+
| lca          | Assigns contig taxonomy based on LCA of all ORF assignments; |
|              | this will be a more stringent and general assignment than    |
|              | ``lca-majority``                                             |
+--------------+--------------------------------------------------------------+
| majority     | Assigns contig taxonomy to tree tip with highest count or    |
|              | tip with highest maximum bitscore                            |
+--------------+--------------------------------------------------------------+

**Default: lca-majority**

::

    annotation:
        references:
            refseq:
                aggregation_method: lca-majority


Majority Threshold
``````````````````

Constitutes a majority fraction for a given tree node within 'lca-majority'
aggregation method.

**Default: 0.51**

::

    annotation:
        references:
            refseq:
                majority_threshold: 0.51


Minimum Alignment Length
````````````````````````

Minimum allowable local alignment length.

**Default: 60**

::

    annotation:
        references:
            refseq:
                min_length: 60


Summary Maximum e-value
```````````````````````

Maximum allowable e-value of high-scoring pair when parsing hits.

**Default: 0.000001**

::

    annotation:
        references:
            refseq:
                max_evalue: 0.000001



Maximum ORF High-Scoring Pairs
``````````````````````````````

The maximum number of hits to consider when summarizing ORFs; can
drastically alter ORF LCA assignments if too high without further limits.

**Default: 10**

::

    annotation:
        references:
            refseq:
                max_hits: 10


Summary Top Fraction of sequences
`````````````````````````````````

Measured from the maximum observed bitscore per query sequence and is used to
further filter after annotation has been completed.

**Default: 0.5**

::

    annotation:
        references:
            refseq:
                top_fraction: 0.5


Minimum Passing Bitscore
````````````````````````

Used to filter hits when parsing after the local alignment step has completed.
A value of zero will disable this filter.

**Default: 0**

::

    annotation:
        references:
            refseq:
                min_bitscore: 500


Example Annotation Section
--------------------------

::

    annotation:
        ## ORFs
        # https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
        translation_table: 11
        # when counting reads aligning to ORFs, require at least this many bp
        # overlapping the ORF
        minimum_overlap: 1
        # count primary read alignments only when getting counts across ORFs
        primary_only: false
        # allow counting of multi-mapped reads when getting counts across ORFs
        multi_mapping: true
        maxbin_max_iteration: 50
        maxbin_min_contig_length: 500
        maxbin_prob_threshold: 0.9
        references:
            eggnog:
                # non-tree based reference requires namemap database and fasta
                namemap: /pic/projects/mint/atlas_databases/eggnog.db
                dmnd: /pic/projects/mint/atlas_databases/eggnog.dmnd
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
                # tree based reference requires namemap database, tree, and diamond formatted fasta
                namemap: /pic/projects/mint/atlas_databases/refseq.db
                tree: /pic/projects/mint/atlas_databases/refseq.tree
                dmnd: /pic/projects/mint/atlas_databases/refseq.dmnd
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
            enzyme:
                namemap: /pic/projects/mint/atlas_databases/enzyme.db
                dmnd: /pic/projects/mint/atlas_databases/enzyme.dmnd
                chunk_size: 500000
                # 'fast' or 'sensitive'
                run_mode: fast
                top_seqs: 2
                index_chunks: 1
                summary_method: majority
            cazy:
                namemap: /pic/projects/mint/atlas_databases/cazy.db
                dmnd: /pic/projects/mint/atlas_databases/cazy.dmnd
                chunk_size: 500000
                # 'fast' or 'sensitive'
                run_mode: fast
                top_seqs: 2
                index_chunks: 1
                summary_method: majority
            cog:
                namemap: /pic/projects/mint/atlas_databases/cog.db
                dmnd: /pic/projects/mint/atlas_databases/cog.dmnd
                chunk_size: 500000
                # 'fast' or 'sensitive'
                run_mode: fast
                top_seqs: 2
                index_chunks: 1
                summary_method: majority
