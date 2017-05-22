Preprocessing of Reads
======================


Adapter Trimming
----------------

FASTA file paths for adapter sequences to be trimmed from the sequence ends.

We provide the adapter reference FASTA included in `bbmap`.

::

    preprocess_adapters: /database_dir/adapters.fa


Quality Trimming
----------------

Trim regions with an average quality below this threshold. Higher is more
stringent.

::

    preprocess_minimum_base_quality: 10


Adapter Trimming at Read Tips
-----------------------------

Allow shorter kmer matches down to `mink` at the read ends. 0 disables.

::

    preprocess_adapter_min_k: 8


Allowable Mismatches in Adapter Hits
------------------------------------

Maximum number of substitutions between the target adapter kmer and the query
sequence kmer. Lower is more stringent.

::

    preprocess_allowable_kmer_mismatches: 1


Contaminant Kmer Length
-----------------------

Kmer length used for finding contaminants. Contaminant matches shorter than
this length will not be found.

::

    preprocess_reference_kmer_match_length: 27


Read Length Threshold
---------------------

This is applied after quality and adapter trimming have been applied to the
sequence.

::

    preprocess_minimum_passing_read_length: 51


Sequence Complexity Filter
--------------------------

Require this fraction of each nucleotide per sequence to eliminate low
complexity reads.

::

    preprocess_minimum_base_frequency: 0.05


Error Correction
----------------

Optionally perform error correction using ``tadpole.sh`` from BBTools::

    perform_error_correction: true


Contamination Parameters
------------------------

Contamination reference sequences in the form of nucleotide FASTA files can be
provided and filtered from the reads using the following parameters::

    contaminant_references:
        # this key ('rRNA') is required
        rRNA: /database_dir/silva_rfam_all_rRNAs.fa
        # any number of arbitrary name:file path references may be added
        phiX: /database_dir/phiX174_virus.fa
    # Don't look for indels longer than this
    contaminant_max_indel: 20
    # Fraction of max alignment score required to keep a site
    contaminant_min_ratio: 0.65
    # mapping kmer length; range 8-15; longer is faster but uses more memory; shorter is more sensitive
    contaminant_kmer_length: 12
    # Minimum number of seed hits required for candidate sites
    contaminant_minimum_hits: 1
    # Set behavior on ambiguously-mapped reads (with multiple top-scoring mapping locations):
    #   best    (use the first best site)
    #   toss    (consider unmapped)
    #   random  (select one top-scoring site randomly)
    #   all     (retain all top-scoring sites)
    contaminant_ambiguous: best


Contaminant References
``````````````````````

As shown in the above example, if provided, reads will be removed from the
FASTQ prior to assembly if they align to these references. If 'rRNA' is
defined, it will be added back to metagenomes prior to assembly.

Additional references can be added arbitrarily, such that::

    contaminant_references:
        rRNA: /refs/rrna.fasta
        human: /refs/human.fasta
        cat: /refs/cat.fasta


Normalization Parameters
````````````````````````

To improve assemblies, coverage is normalized across kmers to a target depth
and can be set using::

    # kmer length over which we calculated coverage
    normalization_kmer_length: 21
    # the normalized target coverage across kmers
    normalization_target_depth: 100
    # reads must have at least this many kmers over min depth to be retained
    normalization_minimum_kmers: 8
