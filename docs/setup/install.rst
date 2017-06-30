Install
=======

All dependencies are installed via conda_ using the bioconda_ channel.
The workflow and some dependencies require Python 3.5::

    conda install python=3.5 snakemake

And install ATLAS::

    pip install pnnl-atlas

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
using `snakemake --use-conda` in the background. No dependencies aside from
python 3.5, snakemake, and pnnl-atlas are required to be installed prior to
execution.

.. _bioconda: https://github.com/bioconda/bioconda-recipes
.. _conda: https://www.continuum.io/downloads
