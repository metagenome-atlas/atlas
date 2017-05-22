Assembly Parameters
===================

Example assembly section::

    assembly:
        assembler: megahit
        memory: 0.99
        minimum_count: 2
        kmer_min: 21
        kmer_max: 121
        kmer_step: 20
        merge_level: 20,0.98
        prune_level: 2
        low_local_ratio: 0.2
        minimum_contig_length: 200
        spades_k: auto
        minc: 5
        minp: 40
        minr: 0
        minl: 200
        trim: 0


Assembler
---------

Currently, the supported assemblers are 'spades' and 'megahit'.

**Default: megahit**

::

    assembly:
        assembler: megahit


``MEGAHIT``: Memory
-------------------

Set the fraction of the machine's total memory if you need to
limit its footprint.

**Default: 0.99**

::

    assembly:
        memory: 0.99


``MEGAHIT``: Minimum Multiplicity
---------------------------------

Set the minimum multiplicity for filtering.

**Default: 2**

::

    assembly:
        minimum_count: 2


``MEGAHIT``: Minimum Kmer Length
--------------------------------

Minimum kmer size (<= 255) and must be odd.

**Default: 21**

::

    assembly:
        kmer_min: 21


``MEGAHIT``: Maximum Kmer Length
--------------------------------

Maximum kmer size (<=255) and must be odd.

**Default: 121**

::

    assembly:
        kmer_max: 121


``MEGAHIT``: Kmer Step
----------------------

Sets the kmer step for ``megahit`` kmer assembly lengths.

**Default: 20**

::

    assembly:
        kmer_step: 20


``MEGAHIT``: Merge Levels
-------------------------

Merge complex bubbles of length <= l*kmer_size and similarity >= s.

**Default: 20,0.98**

::

    assembly:
        merge_level: 20,0.98


``MEGAHIT``: Prune Level
------------------------

Strength of low depth pruning (0-3).

**Default: 2**

::

    assembly:
        prune_level: 2


``MEGAHIT``: Low Local Coverage
-------------------------------

Ratio threshold to define low local coverage contigs.

**Default: 0.2**

::

    assembly:
        low_local_ratio: 0.2


``MEGAHIT``: Minimum Contig Length
----------------------------------

Minimum length of contigs to output from the assembler; can be filtered
downstream using ``minl``.

**Default: 200**

::

    assembly:
        minimum_contig_length: 200


``SPAdes``: Kmer Sizes
----------------------

Comma-separated list of k-mer sizes (must be odd and less than 128).

**Default: auto**

::

    assembly:
        spades_k: auto


Contig Average Coverage Threshold
---------------------------------

Discard contigs with low read support after mapping quality filtered reads
back to contig sequences. Contigs with a lower average coverage than ``minc``
will be removed.

**Default: 5**

::

    assembly:
        minc: 5


Contig Percent Coverage Bases
-----------------------------

Discard contigs with a low fraction of reads mapping back along the length of
the contig.

**Default: 40**

::

    assembly:
        minp: 40


Contig Read Mapping Filter
--------------------------

Require at least this many reads mapped to a contig and discard contigs with
fewer mapped reads.

**Default: 0**

::

    assembly:
        minr: 0


Contig Length Filter
--------------------

Post-assembly contig length filter.

**Default: 1**

::

    assembly:
        minl: 200


Contig Trimming
---------------

Trim the first and last number of bases of each sequence.

**Default: 0**

::

    assembly:
        trim: 0
