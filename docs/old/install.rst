Install
=======

All dependencies are installed via conda_ using the bioconda_ channel.
The workflow and some dependencies require Python 3::

    conda install python>=3.6 snakemake>=5.1.3


The intended usage requires ``conda``.

Using Python 3, install ``atlas``::

    pip install -U metagenome-atlas


Following install, ``atlas`` should be executable::

    $ atlas -h
    Usage: atlas [OPTIONS] COMMAND [ARGS]...

      ATLAS

    Options:
      --version   Show the version and exit.
      -h, --help  Show this message and exit.

    Commands:
      assemble      assembly workflow
      download      download reference files
      make-config   prepopulate a configuration file with samples and defaults
      ...


Execution Environment
---------------------

As ATLAS executes rules to generate output files, an environment is created
using ``snakemake --use-conda`` in the background. No dependencies aside from
Python 3, snakemake, and metagenome-atlas are required to be installed prior to
execution.

For more information related to bioconda, see:
https://bioconda.github.io/

For more information about Snakemake, see:
https://snakemake.readthedocs.io


Execution Paradigm
------------------

The Snakemake environment will be created each time the output directory is
altered. For a given experiment, all samples should be defined in a single
configuration file, so only one environment needs to be created per experiment.

.. _bioconda: https://github.com/bioconda/bioconda-recipes
.. _conda: https://www.continuum.io/downloads
