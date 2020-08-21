


.. _setup_docker:

Docker
******

Use docker container
====================

We recommend to use the conda package as it allows deployment on clusters.
However, if you want to directly start using atlas on a small metagenome you can use the docker container::

  docker pull metagenomeatlas/atlas

Go to a directory on your filesystem where you have the fastq files in a subfolder, e.g. in ``reads``
Your present working directory will be mounted on ``/WD`` in the docker container.

The docker container contains all the dependencies and some of the databases in ``/databases`` .
The databases for functional and taxonomic annotation are downloaded while running.
To not loose the databases after exiting the docker we recommend to mount them also on your disk.

Create::

  mkdir -p AtlasDB/GTDB-TK AtlasDB/EggNOGV2

Then run the docker::

  docker run -i -u $(id -u):$(id -g) -v $(pwd):/WD -v $(pwd)/AtlasDB/EggNOGV2/:/databases/EggNOGV2 -v $(pwd)/AtlasDB/GTDB-TK/:/databases/GTDB-TK -t metagenomeatlas/atlas:latest /bin/bash

Inside the docker you can run atlas as folows::

  atlas init -db-dir /databases /WD/reads

This should create a sample.tsv and a config.yaml, whcih you can edit on your system.
Important don't forget to align the memory of your computer with the memory defined in the config file.

after that run::

  atlas run all


Test atlas using docker container
---------------------------------

Testing atlas using the docker container works similarly to the above example.

Create a working directory::

  mkdir DockerTest; cd DockerTest

get the example data::

  git clone https://github.com/metagenome-atlas/example_data.git

Set up the docker as explained  :ref:`here <docker_setup>` and run it::

  mkdir -p AtlasDB/GTDB-TK AtlasDB/EggNOGV2
  docker run -i -u $(id -u):$(id -g) -v $(pwd):/WD -v $(pwd)/AtlasDB/EggNOGV2/:/databases/EggNOGV2 -v $(pwd)/AtlasDB/GTDB-TK/:/databases/GTDB-TK -t metagenomeatlas/atlas:latest /bin/bash

Inside the docker you have access to the ``example_data`` and you can initialize and run atlas test example as above.
