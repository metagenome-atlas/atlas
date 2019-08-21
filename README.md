# ATLAS

[![Version](https://anaconda.org/bioconda/metagenome-atlas/badges/version.svg)](https://anaconda.org/bioconda/metagenome-atlas)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/metagenome-atlas.svg?label=Bioconda )](https://anaconda.org/bioconda/metagenome-atlas)
[![CircleCI](https://circleci.com/gh/metagenome-atlas/atlas/tree/master.svg?style=svg)](https://circleci.com/gh/metagenome-atlas/atlas/tree/master)
[![Documentation Status](https://readthedocs.org/projects/metagenome-atlas/badge/?version=latest)](https://metagenome-atlas.readthedocs.io/en/latest/?badge=latest)
[![follow on twitter](https://img.shields.io/twitter/follow/SilasKieser.svg?style=social&label=Follow)](https://twitter.com/search?f=tweets&q=%40SilasKieser%20%23metagenomeAtlas&src=typd)
[Slack](https://join.slack.com/t/metagenome-atlas/shared_invite/enQtNTEzMDk2NzI4NjI5LWYxMDVhMzNhMzY3ZDBlOTVjOWI1YzMzNjgwMTZkMDQ0MTNjMDUxZDBhMDkzOTdkMDdiYTAwZDRiOWUwMTY0NDU)

# Quick Start

Three commands to start analysing your metagenome data:
```
    conda install -c bioconda -c conda-forge metagenome-atlas
    atlas init --db-dir databases path/to/fastq/files
    atlas run
```
All databases and dependencies are installed on the fly in the directory `db-dir`.
You want to run these three commands on the example_data on the GitHub repo.
If you have more time, then we recommend you configure atlas according to your needs.
  - check the `samples.tsv`
  - edit the `config.yaml`
  - run atlas on any cluster system
For more details see [documentation](https://metagenome-atlas.rtfd.io/).

# Assembly based metagenomics

Atlas is a easy to use metagenomic pipeline

![scheme of workflow](resources/images/ATLAS_scheme.png?raw=true)


## Install the development version from GitHub
Atlas is still under active development; therefore, you may want to install the up to date atlas from GitHub. Atlas should be run on _linux_, the assembly works also on OS X, but unfortunately not the tools used for binning.
[![Snakemake](https://img.shields.io/badge/snakemake-≥5.4.5-brightgreen.svg)](https://snakemake.bitbucket.io)

Create a conda environment with all primary dependencies. All further dependencies are installed on the fly.
```
conda create -n atlas -c bioconda -c conda-forge python=3.6 snakemake pandas bbmap=37.78 click=7 ruamel.yaml biopython
```
Load the environment:
```
source activate atlas
```
copy code from GitHub and install:
```
git clone https://github.com/metagenome-atlas/atlas.git
cd atlas
pip install --editable .
```
Now you should be able to run atlas:
```
atlas init --db-dir databases path/to/fastq/files
atlas run
```


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
