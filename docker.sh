#!/usr/bin/env bash


echo "
Your present working directory is mounted on /WD in the docker container.
The atlas databases are found on the docker in /databases.

you can run atlas as folows:

  atlas init -db-dir /databases path/to/fastq/inside/your/workindirectory

This should create a sample.tsv and a config.yaml.

after that run:

  atlas run all

"

mkdir -p AtlasDB/GTDB-TK AtlasDB/EggNOGV2
docker run -i -u $(id -u):$(id -g) -v $(pwd):/WD -v AtlasDB/EggNOGV2:/databases/EggNOGV2 -v AtlasDB/GTDB-TK:/databases/GTDB-TK  -t metagenomeatlas/atlas:latest /bin/bash


# https://github.com/metagenome-atlas/example_data
# cd example_data
# docker.sh
# ls
