# Metagenome-Atlas

[![Version](https://anaconda.org/bioconda/metagenome-atlas/badges/version.svg)](https://anaconda.org/bioconda/metagenome-atlas)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/metagenome-atlas.svg?label=Bioconda )](https://anaconda.org/bioconda/metagenome-atlas)
[![CircleCI](https://circleci.com/gh/metagenome-atlas/atlas/tree/master.svg?style=svg)](https://circleci.com/gh/metagenome-atlas/atlas/tree/master)
[![Documentation Status](https://readthedocs.org/projects/metagenome-atlas/badge/?version=latest)](https://metagenome-atlas.readthedocs.io/en/latest/?badge=latest)
[![follow on twitter](https://img.shields.io/twitter/follow/SilasKieser.svg?style=social&label=Follow)](https://twitter.com/search?f=tweets&q=%40SilasKieser%20%23metagenomeAtlas&src=typd)


# Quick Start

Three commands to start analysing your metagenome data:
```
    conda install -y -c bioconda -c conda-forge metagenome-atlas
    atlas init --db-dir databases path/to/fastq/files
    atlas run all
```
All databases and dependencies are installed on the fly in the directory `databases`.

You want to run these three commands on the [example data](https://github.com/metagenome-atlas/atlas/exmple_data).
If you have more time, then we recommend you configure atlas according to your needs.
  - check the `samples.tsv`
  - edit the `config.yaml`
  - run atlas on any [cluster system](https://metagenome-atlas.readthedocs.io/en/latest/usage/cluster.html)
For more details see [documentation](https://metagenome-atlas.rtfd.io/).

# Assembly based metagenomics

Atlas is a easy to use metagenomic pipeline

![scheme of workflow](resources/images/ATLAS_scheme.png?raw=true)

# Setup
Atlas should be run on a _linux_ sytem, with enough memory (min ~50GB but assembly usually requires 250GB).
The only dependency is the _conda package manager_, which can easy be installed with [anaconda](http://anaconda.org/).
We recommend you to create a conda environment for atlas to avoid any conflicts of versions.

```
    conda create -y -n atlasenv
    source activate atlasenv
    conda install -y -c bioconda -c conda-forge metagenome-atlas
```

And you can run atlas. All other dependencies are installed in specific environments during the run of the pipeline.

For local execution we have also a [docker container](https://metagenome-atlas.readthedocs.io/en/latest/usage/getting_started.html#c-use-docker-container)

# License

BSD-3.

# Cite

We have a BioRxiv preprint please cite:

ATLAS: a Snakemake workflow for assembly, annotation, and genomic binning of metagenome sequence data
Silas Kieser, Joseph Brown, Evgeny M Zdobnov, Mirko Trajkovski, Lee Ann McCue
bioRxiv 737528; doi: https://doi.org/10.1101/737528

# Disclaimer

This material was prepared as an account of work sponsored by an agency of the
United States Government.  Neither the United States Government nor the United
States Department of Energy, nor Battelle, nor any of their employees, nor any
jurisdiction or organization that has cooperated in the development of these
materials, makes any warranty, express or implied, or assumes any legal
liability or responsibility for the accuracy, completeness, or usefulness or
any information, apparatus, product, software, or process disclosed, or
represents that its use would not infringe privately owned rights.

Reference herein to any specific commercial product, process, or service by
trade name, trademark, manufacturer, or otherwise does not necessarily
constitute or imply its endorsement, recommendation, or favoring by the United
States Government or any agency thereof, or Battelle Memorial Institute. The
views and opinions of authors expressed herein do not necessarily state or
reflect those of the United States Government or any agency thereof.

PACIFIC NORTHWEST NATIONAL LABORATORY operated by BATTELLE for the UNITED
STATES DEPARTMENT OF ENERGY under Contract DE-AC05-76RL01830
