
.. _longreads:

Long reads
==========

Limitation: Hybrid assembly of long and short reads is supported with spades and metaSpades.
However metaSpades needs a paired-end short-read library.

The path of the (preprocessed) long reads should be added manually to the
the sample table under a new column heading  'longreads'.

In addition the type of the long reads should be defined in the config file:
``longread_type`` one of ["pacbio", "nanopore", "sanger", "trusted-contigs", "untrusted-contigs"]
