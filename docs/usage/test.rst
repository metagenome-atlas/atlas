
.. _example_data:

Test atlas with example data
----------------------------

If you want to test atlas on a small example data here is a two sample, three genome minimal metagenome dataset,
to test atlas. Even when atlas will run faster on the test data,
it will anyway download all the databases and requirements, for the a complete run,
which can take a certain amount of time and especially disk space (>100Gb).

The database dir of the test run should be the same as for the later atlas executions.

The example data can be downloaded as following::

  wget https://zenodo.org/record/3992790/files/test_reads.tar.gz
  tar -xzf test_reads.tar.gz

We initialize a atlas working directory ``testrun`` using the test reads.
The test samples don't require a lot of threads (set to 4), they do require However some memory (~60GB).::

  atlas init --db-dir databases --working-dir testrun --threads 4 test_reads

After the set up you can run::

  atlas run all --working-dir testrun
