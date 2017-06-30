Defining Samples
================

Annotation
----------

Samples are defined under "samples" and have a unique name that does not
contain spaces or underscores (dashes are accepted). To only annotate, a
sample needs ``fasta`` defined::

    samples:
        sample-1:
            fasta: /project/bins/sample-1_contigs.fasta
        sample-2:
            fasta: /project/bins/sample-2_contigs.fasta

All other annotation parameters can be defined following samples. See
:ref:`annotation`.


Annotation and Quantification
-----------------------------

To get counts, we need to define FASTQ file paths that will be mapped back
to the FASTA regions.

In addition to specifying FASTQ paths for each sample, the configuration will
also need to contain::

    quantification: true


Interleaved
```````````

Reads are always assumed to be paired-end, so we only need to specify
``fastq``::

    samples:
        sample-1:
            fasta: /project/bins/sample-1_contigs.fasta
            fastq: /project/data/sample-1_pe.fastq
        sample-2:
            fasta: /project/bins/sample-2_contigs.fasta
            fastq: /project/data/sample-2_pe.fastq

Paired-end
``````````

In this case, we create a list using YAML_ syntax for both R1 and R2 indexes::

    samples:
        sample-1:
            fasta: /project/bins/sample-1_contigs.fasta
            fastq:
                - /project/data/sample-1_R1.fastq
                - /project/data/sample-1_R2.fastq
        sample-2:
            fasta: /project/bins/sample-2_contigs.fasta
            fastq:
                - /project/data/sample-2_R1.fastq
                - /project/data/sample-2_R2.fastq


Single-end
``````````

As data are assumed to be paired-end, we need to add ``paired: false``::

    samples:
        sample-1:
            fasta: /project/bins/sample-1_contigs.fasta
            fastq: /project/data/sample-1_se.fastq
            paired: false
        sample-2:
            fasta: /project/bins/sample-2_contigs.fasta
            fastq: /project/data/sample-2_se.fastq
            paired: false


Example
-------

A complete example for annotation and quantification for samples with
paired-end reads in separate and in interleaved FASTQs::


    samples:
        sample-1:
            fasta: /project/bins/sample-1_contigs.fasta
            fastq: /project/data/sample-1_pe.fastq
        sample-2:
            fasta: /project/bins/sample-2_contigs.fasta
            fastq:
                - /project/data/sample-2_R1.fastq
                - /project/data/sample-2_R2.fastq

    quantification: true

    tmpdir: /scratch
    threads: 24
    refseq_namemap: /pic/projects/mint/atlas_databases/refseq.db
    refseq_tree: /pic/projects/mint/atlas_databases/refseq.tree
    diamond_db: /pic/projects/mint/atlas_databases/refseq.dmnd
    # 'fast' or 'sensitive'
    diamond_run_mode: fast
    # setting top_seqs to 5 will report all alignments whose score is
    # at most 5% lower than the top alignment score for a query
    diamond_top_seqs: 2
    # maximum e-value to report alignments
    diamond_e_value: "0.000001"
    # minimum identity % to report an alignment
    diamond_min_identity: 50
    # minimum query cover % to report an alignment
    diamond_query_coverage: 60
    # gap open penalty
    diamond_gap_open: 11
    # gap extension penalty
    diamond_gap_extend: 1
    # Block size in billions of sequence letters to be processed at a time.
    # This is the main parameter for controlling DIAMOND's memory usage.
    # Bigger numbers will increase the use of memory and temporary disk space,
    # but also improve performance. The program can be expected to roughly use
    # six times this number of memory (in GB).
    diamond_block_size: 6
    # The number of chunks for processing the seed index (default=4). This
    # option can be additionally used to tune the performance. It is
    # recommended to set this to 1 on a high memory server, which will
    # increase performance and memory usage, but not the usage of temporary
    # disk space.
    diamond_index_chunks: 1
    # 'lca', 'majority', or 'best'; summary method for annotating ORFs; when
    # using LCA, it's recommended that one limits the number of hits using a
    # low top_fraction
    summary_method: lca
    # 'lca', 'lca-majority', or 'majority'; summary method for aggregating ORF
    # taxonomic assignments to contig level assignment; 'lca' will result in
    # most stringent, least specific assignments
    aggregation_method: lca-majority
    # constitutes a majority fraction at tree node for 'lca-majority' ORF
    # aggregation method
    majority_threshold: 0.51


.. _YAML: http://www.yaml.org/
