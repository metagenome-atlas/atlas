Annotation
==========

Protocol defaults are reflected in the examples.

Translation Table
-----------------

By default, translation table 11 is used to find open reading frames among
passing contig sequences. Other codes are available at
https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi.

::

    translation_table: 11


Minimum Overlap
---------------

When counting reads overlapping coding sequence, require this much read
overlap.

**Default: 1**

::

    minimum_region_overlap: 1


Restricting Read Counts
-----------------------

Counts can be restricted to primary alignments only using ``primary_only`` in
addition to being able to control behavior associated with multi-mapped reads.
As the alignment stage does allow reads to align ``maximum_counted_map_sites``
per sequence, one may want to later restrict their counts using
``count_multi_mapped_reads``.

::

    primary_only: false
    maximum_counted_map_sites: 10
    count_multi_mapped_reads: false


Genome Binning Options
----------------------

Binning can be skipped entirely by setting ``perform_genome_binning`` to
``false``. If binning is performed, the user can set the following ``maxbin2``
options::

    perform_genome_binning: true
    maxbin_max_iteration: 50
    maxbin_min_contig_length: 500
    maxbin_prob_threshold: 0.9


Functional Annotation of ORFs
-----------------------------

Functional annotation is performed using Prokka. Contigs will be renamed to
sample name + a digit, incrementally, such that contig 1 for sample 'example-id'
is 'example-id_1'. ORFs among a sample are named by Prokka similarly though
they are padded by zeroes (example-id_00001). Contig IDs and ORFs IDs are
mapped back to one another using the final output table where each row
represents an ORF and its assignments.


Taxonomy Annotation of ORFs and Contigs
---------------------------------------

RefSeq version 78 is used for mapping ORFs to products which are then
summarized using NCBI's taxonomy tree. Each ORF is assigned a taxonomy based
on user preference using ``summary_method`` and ``aggregation_method``.
Within the configuration file, a user must define the locations of the RefSeq
files::

    refseq_namemap: /database_dir/refseq.db
    refseq_tree: /database_dir/refseq.tree
    diamond_db: /database_dir/refseq.dmnd


Local Alignment Options for ``blastp`` Search
---------------------------------------------

Within each reference database, the user has the flexibility to optimize
performance across their compute environment and control the number of
alignment hits in various ways.


Run Mode
````````

DIAMOND alignment mode. Either 'fast' of 'more-sensitive'::

    diamond_run_mode: fast


Top Percent of Sequences
````````````````````````

Applies to reported local alignments and will affect the number of hits used
when applying ORF summary methods ``lca`` and ``majority``. Setting
``diamond_top_seqs`` to 5 will report all alignments whose score is at most 5%
lower than the top alignment score for a query.

Set using::

    diamond_top_seqs: 2


Alignment e-value
`````````````````

Maximum e-value to report alignments::

    diamond_e_value: 0.000001


Minimum Identity
````````````````

Filters DIAMOND hits based on minimum matching identity percentage::

    diamond_min_identity: 50


Query Coverage
``````````````

Require this much of the query sequence to be matched above
``diamond_min_identity``::

::

    diamond_query_coverage: 60


Gap Open Penalty
````````````````

A lower gap open penalty may allow more, lower identity hits.

::

    diamond_gap_open: 11


Gap Extend Penalty
``````````````````

A higher extend penalty will reduce allowable indel lengths in matches.

::

    diamond_gap_extend: 1


Block Size
``````````

Block size in billions of sequence letters to be processed at a time.
This is the main parameter for controlling DIAMOND's memory usage.
Bigger numbers will increase the use of memory and temporary disk space,
but also improve performance. The program can be expected to roughly use
six times this number of memory (in GB).

::

    diamond_block_size: 1


Index Chunks
````````````

The number of chunks for processing the seed index. This option can be
additionally used to tune the performance. It is recommended to set this to 1
on a high memory server, which will increase performance and memory usage, but
not the usage of temporary disk space.

::

    diamond_index_chunks: 4


Summary Method
``````````````

This is the summary method for annotating open reading frames. 'lca' performs
an LCA on the hits which can be limited using ``diamond_top_seqs``. Other
options are 'majority' which takes the majority target hit after filtering
alignments and 'best' which simply chooses the top hit.

::

    summary_method: lca


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

::

    aggregation_method: lca-majority


Majority Threshold
``````````````````

Constitutes a majority fraction for a given tree node within 'lca-majority'
aggregation method.

::

    majority_threshold: 0.51
