Test atlas
==========

If you want to test atlas on a small example data here is a two sample, three genome minimal metagenome dataset,
to test atlas. Even when atlas will run faster on the test data,
it will anyway download all the databases and requirements, for the a complete run,
which can take a certain amount of time and especially disk space (>100Gb).

The database dir of the test run should be the same as for the later atlas executions.

The example data can be downloaded as following::

  git clone https://github.com/metagenome-atlas/example_data.git

We initialize a atlas working directory ``testrun`` using the test reads.
The reads are paired end stored in one file therefore we set the ``--interleaved-fastq`` flag.
The test samples don't require a lot of threads (set to 2), they do require However some memory (~60GB).::

  atlas init --db-dir databases --working-dir testrun --interleaved-fastq --threads=2 example_data/reads/test

After the set up you can run::

  atlas run --working-dir testrun


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
