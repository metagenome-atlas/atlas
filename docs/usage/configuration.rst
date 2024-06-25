
..
_configuration:

Configure Atlas
***************

..
_contaminants:

Remove reads from Host
======================

One of the most important steps in the Quality control is to remove reads from the host's genome.
You can add any number of genomes to be removed.

We recommend using genomes where repetitive sequences are masked.
See here for more details `human genome <http://seqanswers.com/forums/archive/index.php/t-42552.html>`_.


Co-abundance Binning
====================

.. _cobinning:

While binning each sample individually is faster, using co-abundance for binning is recommended.
Quantifying the coverage of contigs across multiple samples provides valuable insights about contig co-variation.

There are two primary strategies for co-abundance binning:

1. **Cross mapping:** Map the reads from multiple samples to each sample's contigs.
2. **Co-binning:** Concatenate contigs from multiple samples and map all the reads to these combined contigs.

`final_binner: metabat2` is used for cross-mapping, while `vamb` or `SemiBin` is used for co-binning.

The samples to be binned together are specified using the `BinGroup` in the `sample.tsv` file.
The size of the BinGroup should be selected based on the binner and the co-binning strategy in use.

Cross-mapping complexity scales quadratically with the size of the BinGroup since each sample's reads are mapped to each other.
This might yield better results for complex metagenomes, although no definitive benchmark is known.
On the other hand, co-binning is more efficient, as it maps a sample's reads only once to a potentially large assembly.

Default Behavior
----------------

Starting with version 2.18, Atlas places every sample in a single BinGroup and defaults to `vamb` as the binner unless there are very few samples.
For fewer than 8 samples, `metabat` is the default binner.

.. note::
    This represents a  departure from previous versions, where each sample had its own BinGroup.
    Running `vamb` in those versions would consider all samples, regardless of their BinGroup.
    This change might cause errors if using a `sample.tsv` file from an older Atlas version.
    Typically, you can resolve this by assigning a unique BinGroup to each sample.

The mapping threshold has been adjusted to 95% identity (single sample binning is 97%) to allow reads from different strains — 
but not other species — to map to contigs from a different sample.

If you're co-binning more than 150-200 samples or cross-mapping more than 50 samples, Atlas will issue a warning regarding excessive samples in a BinGroup.
Although VAMB's official publication suggests it can handle up to 1000 samples, this demands substantial resources.

Therefore, splitting your samples into multiple BinGroups is recommended.
Ideally, related samples, or those where the same species are anticipated, should belong to the same BinGroup.

Single-sample Binning
---------------------

To employ single-sample binning, simply assign each sample to its own BinGroup and select `metabat` or `DASTool` as the `final_binner`.

Although it's not recommended, it's feasible to use `DASTool` and feed it inputs from `metabat` and other co-abundance-based binners.

Add the following lines to your `config.yaml`:


.. code-block:: yaml

   final_binner: DASTool

   binner: 
     - metabat
     - maxbin
     - vamb



.. _longreads:

Long reads
==========

Limitation: Hybrid assembly of long and short reads is supported with spades and metaSpades.
However, metaSpades needs a paired-end short-read library.

The path of the (preprocessed) long reads should be added manually to the
sample table under a new column heading  'longreads'.

In addition, the type of the long reads should be defined in the config file:
``longread_type`` one of ["pacbio", "nanopore", "sanger", "trusted-contigs", "untrusted-contigs"]


Example config file
===================


..include:: ../../config/template_config.yaml
  :code:




Detailed configuration
======================

..
toctree::
    :maxdepth: 1

    ../advanced/qc
    ../advanced/assembly
