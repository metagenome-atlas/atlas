#!/usr/bin/env bash

curdir=$(pwd)

databases=./databases # path to database dir where you can store some 100gb databases

# install atlas
git clone https://github.com/metagenome-atlas/atlas.git
cd atlas
conda create f atlasenv.yml
conda activate atlasenv
pip install --editable .


# set up cluster integration
mkdir -p ~/.config/snakemake
cd ~/.config/snakemake
cookiecutter https://github.com/metagenome-atlas/clusterprofile.git

cd $curdir
#test run
git clone https://github.com/metagenome-atlas/example_data.git
atlas init --db-dir $databases --working-dir testrun example_data/reads/test
atlas run all --working-dir testrun --profile cluster
