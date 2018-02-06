Pre-Assambly-processing
=======================

Normalization Parameters
------------------------

To improve assembly time and often assemblies themselves, coverage is
normalized across kmers to a target depth and can be set using::

    # kmer length over which we calculated coverage
    normalization_kmer_length: 21
    # the normalized target coverage across kmers
    normalization_target_depth: 100
    # reads must have at least this many kmers over min depth to be retained
    normalization_minimum_kmers: 8



Error Correction
----------------

Optionally perform error correction using ``tadpole.sh`` from BBTools::

    perform_error_correction: true



Assembly Parameters
===================


Assembler
---------

Currently, the supported assemblers are 'spades' and 'megahit' with the
default setting of::

    assembler: megahit

Both assemblers have settings that can be altered in the configuration::

    # minimum multiplicity for filtering (k_min+1)-mers
    megahit_min_count: 2
    # minimum kmer size (<= 255), must be odd number
    megahit_k_min: 21
    # maximum kmer size (<= 255), must be odd number
    megahit_k_max: 121
    # increment of kmer size of each iteration (<= 28), must be even number
    megahit_k_step: 20
    # merge complex bubbles of length <= l*kmer_size and similarity >= s
    megahit_merge_level: 20,0.98
    # strength of low depth pruning (0-3)
    megahit_prune_level: 2
    # ratio threshold to define low local coverage contigs
    megahit_low_local_ratio: 0.2
    # minimum length of contigs (after contig trimming)
    minimum_contig_length: 200
    # comma-separated list of k-mer sizes (must be odd and less than 128)
    spades_k: auto


Contig Filtering
----------------

After assembly, contigs can be filtered based on several metrics::

    # Discard contigs with lower average coverage.
    minimum_average_coverage: 5
    # Discard contigs with a lower percent covered bases.
    minimum_percent_covered_bases: 40
    # Discard contigs with fewer mapped reads.
    minimum_mapped_reads: 0
    # Trim the first and last X bases of each sequence.
    contig_trim_bp: 0
