ATLAS
=====


Quick Start
-----------

Three commands to start analyzing your metagenome data::
    conda install -c bioconda -c conda-forge metagenome-atlas
    atlas init --db-dir databases path/to/fastq/files
    atlas run

All databases and dependencies are installed on the fly in the directory ``--db-dir``.
You want to run this three commands on the example_data on the GitHub repo.
If you have more time, then  we recommend you to configure atlas according to your needs.
  - check the ``samples.tsv``
  - edit the ``config.yaml``
  - run atlas on a cluster system
For more details see documentation.

|scheme|


Documentation
=============

.. toctree::
    :maxdepth: 2
    :caption: Setup

    usage/getting_started
    usage/threads

.. toctree::
    :maxdepth: 2
    :caption: Outputs

    pipeline

.. toctree::
    :maxdepth: 1
    :caption: Advanced configuration

    advanced/qc
    advanced/assembly


Indices and tables
==================

.. |scheme| image:: ../resources/images/ATLAS_scheme.png
  :alt: Atlas is a workflow for assembly and binning of metagenomic reads


* :ref:``search``
