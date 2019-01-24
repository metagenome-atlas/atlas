ATLAS
=====


Quick Start
-----------

Three commands to start analyzing your metagenome data::
    conda install -c bioconda -c conda-forge metagenome-atlas #may not yet be ready see below
    atlas init --db-dir databases path/to/fastq/files
    atlas run

All databases and dependencies are installed on the fly in the directory ``db-dir``.
You want to run this three commands on the example_data on the GitHub repo.
If you have more time, then  we recommend you to configure atlas according to your needs.
  - check the ``samples.tsv``
  - edit the ``config.yaml``
  - run atlas on a cluster system
For more details see documentation.




.. toctree::
    :maxdepth: 2
    :caption: Install

    getting_started
    threads


.. toctree::
    :maxdepth: 2
    :caption: Assembly Protocol

    assembly/samples
    assembly/threads
    assembly/preprocessing
    assembly/assembly
    assembly/annotation
    assembly/advanced
    assembly/example_configuration
    assembly/output


.. toctree::
    :maxdepth: 2
    :caption: Annotation Protocol

    annotation/samples


Indices and tables
==================

* :ref:``search``
