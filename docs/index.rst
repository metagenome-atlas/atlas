.. image:: https://img.shields.io/conda/dn/bioconda/metagenome-atlas.svg?label=Bioconda
    :target: https://bioconda.github.io/recipes/metagenome-atlas/README.html

.. .. image:: https://img.shields.io/pypi/pyversions/metagenome-atlas.svg
    :target: https://www.python.org

.. image:: https://img.shields.io/pypi/v/metagenome-atlas.svg
    :target: https://pypi.python.org/pypi/metagenome-atlas

.. .. image:: https://quay.io/repository/snakemake/snakemake/status
       :target: https://quay.io/repository/snakemake/snakemake

.. .. image:: https://circleci.com/bb/metagenome-atlas/atlas/tree/master.svg?style=shield
    :target: https://circleci.com/bb/metagenome-atlas/atlas/tree/master

.. image:: https://img.shields.io/badge/stack-overflow-orange.svg
    :target: http://stackoverflow.com/questions/tagged/metagenome-atlas

.. image:: https://img.shields.io/twitter/follow/SilasKieser.svg?style=social&label=Follow
    :target: https://twitter.com/search?f=tweets&q=%40SilasKieser%20%23metagenomeAtlas&src=typd


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
If you have more time, then  we recommend you to configure atlas according to your needs:

  - check the samples.tsv
  - edit the :ref:`config.yaml <configuration>`
  - run atlas on a :ref:`cluster system <execution_system>`

|scheme|


Documentation
=============

.. toctree::
    :maxdepth: 1
    :caption: Setup

    usage/getting_started
    usage/threads

.. toctree::
    :maxdepth: 2
    :caption: The Pipeline

    details/pipeline
    details/configuration





.. |scheme| image:: ../resources/images/ATLAS_scheme.png
  :alt: Atlas is a workflow for assembly and binning of metagenomic reads
