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

::

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



annotation:
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
        enzyme:
            namemap: /pic/projects/mint/atlas_databases/functional/enzyme/enzyme.db
            fasta: /pic/projects/mint/atlas_databases/functional/enzyme/enzyme.fasta
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
