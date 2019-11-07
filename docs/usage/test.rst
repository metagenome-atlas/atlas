Test atlas
==========

If you want to test atlas on a small example data here is a two sample, three genome minimal metagenome dataset,
to test atlas. Even when atlas will run faster on the test data,
it will anyway download all the databases and requirements, for the a complete run,
which can take a certain amount of time and especially disk space (>100Gb).

The database dir of the test run should be the same as for the later atlas executions.

The example data can be downloaded as following::

  git clone https://github.com/metagenome-atlas/example_data.git

  atlas init --db-dir databases --working-dir testrun example_data/reads/test

  atlas run --working-dir testrun
