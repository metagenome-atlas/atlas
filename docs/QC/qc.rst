Quality control of reads
========================


Adapter Trimming
----------------

FASTA file paths for adapter sequences to be trimmed from the sequence ends.

We provide the adapter reference FASTA included in `bbmap` for various

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


Contamination Parameters
------------------------

Contamination reference sequences in the form of nucleotide FASTA files can be
provided and filtered from the reads using the following parameters.

If 'rRNA' is defined, it will be added back to metagenomes but not to metatranscriptomes.
Additional references can be added arbitrarily, such as::
::

    contaminant_references:
        rRNA: /database_dir/silva_rfam_all_rRNAs.fa
        phiX: /database_dir/phiX174_virus.fa

Don't look for indels longer than this::

    contaminant_max_indel: 20


Fraction of max alignment score required to keep a site::

    contaminant_min_ratio: 0.65
mapping kmer length; range 8-15; longer is faster but uses more memory; shorter is more sensitive::

    contaminant_kmer_length: 12

Minimum number of seed hits required for candidate sites::

    contaminant_minimum_hits: 1

Set behavior on ambiguously-mapped reads (with multiple top-scoring mapping locations):

- best    (use the first best site)
- toss    (consider unmapped, retain in reads for assembly)
- random  (select one top-scoring site randomly)
- all     (retain all top-scoring sites)

::

    contaminant_ambiguous: best


For human reads we recommend to download the genome [here](), where contaminants and low complexity regions were masked.
