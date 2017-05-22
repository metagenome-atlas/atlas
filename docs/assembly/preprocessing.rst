Preprocessing of Reads
======================

All settings will be indented from 'preprocessing' and indentation levels
matter. See :ref:`example-configuration` for details.

A minimal example of the preprocessing section::

    preprocessing:
        adapters: /databases/adapters.fa
        minimum_base_quality: 10
        min_base_frequency: 0.05
        contamination:
            references:
                rRNA: /refs/rrna.fasta
            k: 12
            ambiguous: best
        normalization:
            k: 21
            t: 100


Adapters
--------

FASTA file paths for adapter sequences to be trimmed from the sequence ends.
It is best practice to use full file paths to reference files.

We provide the adapter reference FASTA included in `bbmap`.

::

    preprocessing:
        adapters: /databases/adapters.fa


Quality Trimming
----------------

Trim regions with an average quality below this threshold. Higher is more
stringent.

**Default: 10**

::

    preprocessing:
        minimum_base_quality: 10


Adapter Trimming at Read Tips
-----------------------------

Allow shorter kmer matches down to `mink` at the read ends. 0 disables.

**Default: 8**

::

    preprocessing:
        mink: 8


Allowable Mismatches in Adapter Hits
------------------------------------

Maximum number of substitutions between the target adapter kmer and the query
sequence kmer. Lower is more stringent.

**Default: 1**

::

    preprocessing:
        allowable_kmer_mismatches: 1


Contaminant Kmer Length
-----------------------

Kmer length used for finding contaminants. Contaminant matches shorter than
this length will not be found.

**Default: 27**

::

    preprocessing:
        reference_kmer_match_length: 27


Read Length Threshold
---------------------

This is applied after quality and adapter trimming have been applied to the
sequence.

**Default: 51**

::

    preprocessing:
        minimum_passing_read_length: 51


Sequence Complexity Filter
--------------------------

Require this fraction of each nucleotide per sequence to eliminate low
complexity reads.

**Default: 0.05**

::

    preprocessing:
        min_base_frequency: 0.05


Contamination Parameters
------------------------

Contamination reference sequences in the form of nucleotide FASTA files can be
provided and filtered from the reads using the following parameters. These
still fall within the 'preprocessing' section of the configuration.


Maximum Insertion/Deletion
``````````````````````````

Have `bbsplit.sh` stop searching for possible mappings with indels longer than
this. Lower is faster.

**Default: 20**

::

    preprocessing:
        contamination:
            maxindel: 20



Required Mapped Read Fraction
`````````````````````````````

Of the possible maximum alignment score, force at least this fraction per mapping.

**Default: 0.65**

::

    preprocessing:
        contamination:
            minratio: 0.65


Minimum Seed Hits
`````````````````

Minimum number of seed hits required for candidate sites.

**Default: 1**

::

    preprocessing:
        contamination:
            minhits: 1


Ambiguous Mappings
``````````````````

The method for which we will deal with reads that map to multiple
contamination reference sequences. Possible values include:

+--------+---------------------------------------+
| Value  | Definition                            |
+========+=======================================+
| best   | Use the first best site.              |
+--------+---------------------------------------+
| toss   | Consider the read unmapped.           |
+--------+---------------------------------------+
| random | Select one top-scoring site randomly. |
+--------+---------------------------------------+
| all    | Retain all top-scoring sites.         |
+--------+---------------------------------------+

**Default: best**

::

    preprocessing:
        contamination:
            ambiguous: best


Mapping Kmer Length
```````````````````

Mapping kmer length in the range of 8 to 15. Shorter will be more sensitive
and slower.

**Default: 13**

::

    preprocessing:
        contamination:
            k: 13


Reference Sequences
```````````````````

Reference FASTA files are defined under 'contamination'.


Ribosomal RNA (``rRNA``)
''''''''''''''''''''''''''''

This reference FASTA is required though you can provide an alternate to the
provided rRNA reference.

::

    preprocessing:
        contamination:
            references:
                rRNA: /refs/rrna.fasta


Additional References
'''''''''''''''''''''

Any number of additional contamination reference sequences can be used. The
key is the name that will be integrated into the file name and provide the
path to the file such that:

::

    preprocessing:
        contamination:
            references:
                rRNA: /refs/rrna.fasta
                human: /refs/human.fasta
                cat: /refs/cat.fasta


Normalization Parameters
````````````````````````

To improve assemblies, coverage is normalized across kmers to a target depth.


Kmer Length
'''''''''''

Kmer length over which we calculated coverage.

**Default: 21**

::

    preprocessing:
        normalization:
            k: 21


Target Coverage
'''''''''''''''

The normalized target coverage across kmers.

**Default: 100**

::

    preprocessing:
        normalization:
            t: 100


Minimum Passing Kmers
'''''''''''''''''''''

Reads must have at least this many kmers over the minimum depth to be retained.

**Default: 8**

::

    preprocessing:
        normalization:
            minkmers: 8


Full Example with Advanced Options
----------------------------------

::

    preprocessing:

        adapters: /pic/projects/mint/atlas_databases/adapters.fa
        mink: 8
        minimum_base_quality: 10
        allowable_kmer_mismatches: 1
        reference_kmer_match_length: 31
        minimum_passing_read_length: 51
        min_base_frequency: 0.05

        contamination:
            references:
                rRNA: /pic/projects/mint/atlas_databases/silva_rfam_all_rRNAs.fa
                phiX: /pic/projects/mint/atlas_databases/phiX174_virus.fa
            maxindel: 20
            minratio: 0.65
            k: 12
            minhits: 1
            ambiguous: best

        normalization:
            k: 21
            t: 100
            minkmers: 8
