
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


Co-abundance binning
======================

.. _cobinning:

While it is faster to bin each sample individually, we recommend using co-abundance for binning.
By quantifying the coverage of contigs in multiple samples one can obtain useful information about which contigs co-vary across samples.

There are two strategies to use co-abundance for binning.

1. Cross mapping: mapping the reads of multiple samples to each samples' contigs.
2. Co-binning: Concatenate the contigs of multiple samples and map the all the reads to the combined contigs.

Cross-mapping is done with `final_binner: metabat2`, while Co-binning is done with `vamb`, or `SemiBin`.

Which samples are binned together are defined with the `BinGroup` in the sample.tsv.
The Size of the BinGroup should be chosen depending on the binner, respectively the co-bining strategy.

Cross mapping scales quadratically with the size of the BinGroup, as each samples' reads are mapped to each other.
It might yeld better results for very complicated metagenomes, even though I am not aware of a good benchmark. 
Co-binning is more efficient as it requires a samples' reads to be mapped only to one, even tough large, assembly. 


If you have less than 5 



You need to put multiple samples in the same `BinGroup` in the sample.tsv.

Atlas requires to have at least 5 samples in a group to use co-abundance binning.
It also warns you not to have too many samples, because then it will take too long to run. 

Since atlas. v2.18 the co-binni




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
