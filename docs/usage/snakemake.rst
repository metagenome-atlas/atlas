.. _snakemake:


Snakemake command line arguments
================================

Atlas builds on snakemake. We designed the command line interface in a way that additional snakemake arguments can be added to an atlas run call.

For instance the ```--profile`` used for cluster execution. Other handy snakemake command line arguments include.

 ``--keep-going``, which  allows atlas in the case of a failed job to continue with independent steps.

For a full list of snakemake arguments see the `snakemake doc <https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options>`_.
