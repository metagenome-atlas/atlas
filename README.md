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

# Citation

ATLAS: a Snakemake workflow for assembly, annotation, and genomic binning of metagenome sequence data.
Kieser, S., Brown, J., Zdobnov, E. M., Trajkovski, M. & McCue, L. A. 
BMC Bioinformatics 21, 257 (2020).
doi: [10.1186/s12859-020-03585-4](https://doi.org/10.1186/s12859-020-03585-4)


