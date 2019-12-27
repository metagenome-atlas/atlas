#!/usr/bin/env bash

databases=./databases # path to database dir where you can store some 100gb databases

# install atlas
git clone https://github.com/metagenome-atlas/atlas.git
cd atlas
conda env create -n atlasenv -f atlasenv.yml
conda activate atlasenv
pip install --editable .

atlas --help

#test run
git clone https://github.com/metagenome-atlas/example_data.git
atlas init --db-dir $databases --working-dir testrun --interleaved-fastq example_data/reads/test

# set up cluster integration
mkdir -p ~/.config/snakemake
cookiecutter https://github.com/metagenome-atlas/clusterprofile.git -o ~/.config/snakemake

atlas run all --working-dir testrun --profile cluster --dryrun
atlas run all --working-dir testrun --profile cluster
