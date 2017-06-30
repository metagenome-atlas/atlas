Getting Started
===============

After installing, one needs to download the required databases and
create a sample configuration file.

Databases
---------

To download the databases and their respective metadata databases:

::

    atlas download -o ~/databases

The downloads use approximately 30 GB of disk space.

Configuration File
------------------

To create a configuration file run:

::

    atlas make-config --database-dir ~/databases \
        config.yaml ~/directory_with_fastqs

Sample names and file paths along with default settings will populate
config.yaml. This `YAML <http://www.yaml.org/start.html>`__ file can be
updated with any text editor.

Sample names should be A-Z characters and can be dash ("-") delimited.

Assembly
--------

After editing your configuration file and adjusting any additional
parameters we run assemblies across our samples using:

::

    atlas assemble config.yaml

By default, this will write results into our current working directory
across the total number of CPU cores available.
