Install
=======

All dependencies are installed via conda_ using the bioconda_ channel.
The workflow and some dependencies require Python 3.5.

To install::

    conda install -c bioconda \
        bbmap diamond fastqc megahit prodigal samtools snakemake spades verse


Or as an isolated environment using our `environment.yml` file::

    conda env create -f environment.yml


And load and unload that environment using ``source activate atlas_env``
and ``source deactivate atlas_env``, respectively.

In the future we plan to push ``atlas`` to Bioconda and PyPI, but currently
to install you will need to download or clone ``atlas``.

Then within that environment and from within the 'atlas' folder::

    cd atlas
    python setup.py install

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


.. _bioconda: https://github.com/bioconda/bioconda-recipes
.. _conda: https://www.continuum.io/downloads
