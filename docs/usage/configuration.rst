
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
It might yeld better results for very complicated metagenomes, but I am not aware of a good benchmark. 
Co-binning is more efficient as it requires a samples' reads to be mapped only to one, even tough large, assembly. 

default behavior
`````````````````

Since version 2.18 atlas, puts every sample in one BinGroup and sets's vamb as default binner unless you have very view samples.
If you have less than 8 samples it uses metabat as binner.

*Note: This is a small but breaking change to previous versions where all samples where in their own BinGroup, and running vamb resulted in taking all samples regardless of their BinGroup.
This means you might see errors when using a sample.tsv generated with an older version of atlas. In most cases you simply need to set a unique BinGroup for all samples.*

Also of note is that we adapted the mapping threshold to 95% identity (vs 97%) to allow the mapping of reads from different strains,
 but not other species to map to contigs from different sample.

If you have more than 150-200 samples for co-binning or more than 50 samples for cross-mapping, atlas will print a warning that too many samples are in a BinGroup.
While according to the official publication of vamb, it can be run of up to 1000 samples, this requires very high resources.

We therefore recommend to split your samples into BinGroups. Ideally samples that are more or less related, where you expect the same species should be in the same BinGroup. 

Single-sample Binning
`````````````````````

If you really want to use single-sample binning. You simply need to put each sample in it's own binGroup and use `metabat` or `DASTool` as your final_binner.

It is not recommend, but possible to use DASTool, and give it input from metabat and other binners based on co-abundance. 

You can add the folowing lines to your config.yaml.

```
final_binner: DASTool

binner: 
  - metabat
  - maxbin
  - vamb
```



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
