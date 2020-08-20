.. image:: https://img.shields.io/conda/dn/bioconda/metagenome-atlas.svg?label=Bioconda
    :target: https://bioconda.github.io/recipes/metagenome-atlas/README.html

.. image:: https://img.shields.io/twitter/follow/SilasKieser.svg?style=social&label=Follow
    :target: https://twitter.com/search?f=tweets&q=%40SilasKieser%20%23metagenomeAtlas&src=typd


Metagenome-Atlas
================

Metagenome-atlas is a easy-to-use metagenomic pipeline based on [snakemake](https://snakemake.github.io/). It handles all steps from QC, Assembly, Binning, to Annotation.

|scheme|

.. |scheme| image:: ../../resources/images/ATLAS_scheme.png
  :alt: Atlas is a workflow for assembly and binning of metagenomic reads


You can start using atlas with three commands::

      conda install -c bioconda -c conda-forge metagenome-atlas
      atlas init --db-dir databases path/to/fastq/files
      atlas run




Documentation
=============

.. toctree::
    :maxdepth: 1
    :caption: Getting started

    usage/getting_started
    usage/test
    usage/cluster
    usage/docker


.. toctree::
    :maxdepth: 2
    :caption: Documentation

    details/pipeline
    details/longreads
    details/configuration
