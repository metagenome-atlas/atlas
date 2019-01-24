Assembly Output
===============

Annotations and contigs are found in the top-level of the output. The contigs
represent the passing contigs after filtering has been performed. The
annotations table contains metrics on the contigs, open reading frames,
genomics bins, taxonomy assignments, other functional annotations, and the
coverage count across the genomic regions.


Sequence Quality Control
------------------------

Paired-end and single-end reads with 00 have been quality trimmed and filtered
of adapter sequences using BBDuk2. The minimum length threshold and sequence
complexity filters are also applied to these reads.

Paired-end reads with 01 represent error corrected reads using Tadpole of the
BBTools package.

FASTQ files with 02 represent sequences pulled from 01, or 00 if no error
correction has been performed, that are hits to the specified contamination
reference database defined in the file name.

Sequence reads in 03 represent reads after decontamination where if the
sample is a metagenome, rRNA reads have been added back into the pool. rRNA
reads are not added back for metatranscriptomes.

Reads that have normalized coverage across kmers are in 04.


Assembly
--------

Unfiltered contigs and stats associated with the assembly, coverage across
contigs before and after filtering, and logs associated with these steps.

Prefiltered contigs represents unfiltered contigs immediately following
assembly.

Discarded contigs represents a subset of contigs present in prefiltered
contigs that were omitted from the final contigs due to coverage or length
parameters set by the user.


Sequence Alignment
------------------

Holds temporary alignment files as well as the final alignment file after
Picard MarkDuplicates has been run against the alignments. This file is used
in the quantification step by ``featureCounts``.


Annotation
----------

Holds functional annotation output and raw counts across functional regions.

Feature Counts
``````````````

Per loci read counts and a summary of read mapping for the sample.

Prokka
``````

All Prokka output including the amino acid sequence FASTA and associated
functional annotation. The "plus" file represents a fixed TSV format that
ensures all rows have the same number of columns.

RefSeq
``````

Includes raw BLAST hits, sorted BLAST hits with an added first column which
enables contig to open reading frame mapping, and taxonomy assignments per
loci.


Genomic Bins
------------

Contains the bins determined by MaxBin2 and all checkm output including its
taxonomy assignments and completeness estimates per bin. These metrics are
combined into a single table in the top level of the sample output directory.


Reference Genome (ref)
----------------------

Contains a concatenated genome of contamination references and will be used
on subsequent samples being output to the same results output directory.
