
Quick Start
===========

Three commands to start analyzing your metagenome data::

    conda install -c bioconda -c conda-forge metagenome-atlas
    atlas init --db-dir databases path/to/fastq/files
    atlas run

All databases and dependencies are installed on the fly in the directory `db-dir`.
You want to run this three commands on the example_data on the GitHub repo.
If you have more time, then  we recommend to take care on the following steps.

Step-by-step installation
=========================

1. Create conda environment
---------------------------

You need to install anaconda or miniconda.
We recommend you to create a conda environment::

    conda create -n atlasenv
    conda activate atlasenv

Then install metagenome-atlas::

    conda install -y -c bioconda -c conda-forge metagenome-atlas


2. Download all databases first
-------------------------------

May be you want to make sure that all databases are downloaded correctly. Simply run::

    atlas download --db-dir path/to/databases

To reassure you, most of the databases are md5 checked. The downloads use approximately 30 GB of disk space.

3. Test installation
--------------------

Use our example_data on the GitHub repo. The first time you run atlas, it installs all dependecies.
It needs therefore an internet connection and some time.

Usage
=====

Now let's apply atlas on your data.

atlas init
----------

::

    Usage: atlas init [OPTIONS] PATH_TO_FASTQ

      Write the file CONFIG and complete the sample names and paths for all
      FASTQ files in PATH.

      PATH is traversed recursively and adds any file with '.fastq' or '.fq' in
      the file name with the file name minus extension as the sample ID.

    Options:
      -d, --db-dir PATH               location to store databases (need ~50GB)
                                      [default: /Users/silas/Documents/Debug_atlas
                                      /databases]
      -w, --working-dir PATH          location to run atlas
      --assembler [megahit|spades]    assembler  [default: megahit]
      --data-type [metagenome|metatranscriptome]
                                      sample data type  [default: metagenome]
      --threads INTEGER               number of threads to use per multi-threaded
                                      job
      -h, --help                      Show this message and exit.


This command creates a ``samples.tsv`` and a ``config.yaml`` in the working directory.

Have a look at them with a normal text editor and check if the samples names are inferred correctly.
Samples should be alphanumeric names and cam be dash delimited. Underscores should be fine too.

You can find an example here. ADDLINK

The ``BinGroup`` parameter is used during the genomic binning.
In short: all samples in which you expect the same strain to
be found should belong to the same group,
e.g. all metagenome samples from mice in the same cage. ADDLINK

You should also check the ``config.yaml`` file, especially:

- You may want to change the resources configuration, depending on the system you run atlas. ADDLINK
- You may want to add ad host genomes to be removed. ADDLINK

Details about the parameters can be found here. ADDLINK

atlas run
----------

::

  Usage: atlas run [OPTIONS] [[qc|assembly|genomes|genecatalog|None|all]]
                   [SNAKEMAKE_ARGS]...

    Runs the ATLAS pipline

    By default all steps are executed but a sub-workflow can be specified.
    Needs a config-file and expects to find a sample table in the working-
    directory. Both can be generated with 'atlas init'

    Most snakemake arguments can be appended to the command for more info see
    'snakemake --help'

    For more details, see: https://metagenome-atlas.readthedocs.io

  Options:
    -w, --working-dir PATH  location to run atlas.
    -c, --config-file PATH  config-file generated with 'atlas init'
    -j, --jobs INTEGER      use at most this many jobs in parallel (see cluster
                            submission for mor details).  [default: 8]
    --no-conda              do not use conda environments. good luck!  [default:
                            False]
    -n, --dryrun            Test execution.  [default: False]
    -h, --help              Show this message and exit.


``atlas run`` need to know the working directory with a ``samples.tsv`` inside it.

Take note of the ``--dryrun`` parameter, see here for other handy snakemake arguments. ADDLINK

If you want to run atlas on a cluster system you want to read this section.
