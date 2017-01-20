Getting Started
===============

After installing, one needs to download the required databases and create a
sample configuration file.


Databases
---------

To download the databases and their respective metadata databases::

    atlas download -o /databases

The downloads use approximately 40 GB of disk space.


Configuration File
------------------

To create a simple configuration file (which you can later edit), run::

    atlas make-config --database-dir /databases config.yaml /my-data

In this case, 'my-data' is a directory containing .fastq files similar to::

    $ tree /my-data
    /my-data/
    ├── Sample-1_R1.fastq.gz
    ├── Sample-1_R2.fastq.gz
    ├── Sample-2_R1.fastq.gz
    └── Sample-2_R2.fastq.gz

Paired-end, single-end (currently with some caveats), and interleaved paired-end
FASTQs are supported.

.. important::
    Single-end currently works, but parameters to assemblers still need to be
    altered such that they are not input as interleaved paired-end.


Assembly
--------

After editing your configuration file and adjusting any additional parameters
we run assemblies across our samples using::

    atlas assemble config.yaml

By default, this will write results into our current working directory across
the total number of CPU cores available.
