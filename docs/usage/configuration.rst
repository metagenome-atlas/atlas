
.. _configuration:

Configure Atlas
***************

.. _contaminants:

Remove reads from Host
======================

One of the most important steps in the Quality control is to remove host genome.
You can add any number of genomes to be removed.

We recommend you to use genomes where repetitive sequences are masked.
See here for more details `human genome <http://seqanswers.com/forums/archive/index.php/t-42552.html>`_.


.. _longreads:

Long reads
==========

Limitation: Hybrid assembly of long and short reads is supported with spades and metaSpades.
However metaSpades needs a paired-end short-read library.

The path of the (preprocessed) long reads should be added manually to the
the sample table under a new column heading  'longreads'.

In addition the type of the long reads should be defined in the config file:
``longread_type`` one of ["pacbio", "nanopore", "sanger", "trusted-contigs", "untrusted-contigs"]


Example config file
===================


.. include:: ../../workflow/../config/template_config.yaml
  :code:




Detailed configuration
======================

.. toctree::
    :maxdepth: 1

    ../advanced/qc
    ../advanced/assembly
