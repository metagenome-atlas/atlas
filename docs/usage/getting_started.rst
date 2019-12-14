Install
========

A. Use conda
-------------

You need to install [anaconda](http://anaconda.org/) or miniconda.
We recommend you to create a conda environment, then install metagenome-atlas::

    conda create -y -n atlasenv
    source activate atlasenv
    conda install -y -c bioconda -c conda-forge metagenome-atlas


B. Install the development version from GitHub
-----------------------------------------------
Atlas is still under active development, therefore you may want to install the up to date atlas from GitHub.

get code from GitHub::

  git clone https://github.com/metagenome-atlas/atlas.git
  cd atlas

Create a conda environment with all primary dependencies. All further dependencies are installed on the fly::

  conda env create -f atlasenv.yml
  source activate atlasenv

Install atlas::

  pip install --editable .


Now you should be able to run atlas::

  atlas init --db-dir databases path/to/fastq/files
  atlas run

C. Use docker container
-----------------------

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

  docker run -i -u $(id -u):$(id -g) -v $(pwd):/WD -v AtlasDB/EggNOGV2:/databases/EggNOGV2 -v AtlasDB/GTDB-TK:/databases/GTDB-TK  -t metagenomeatlas/atlas:latest /bin/bash


Inside the docker you can run atlas as folows::

  atlas init -db-dir /databases /WD/reads

This should create a sample.tsv and a config.yaml, whcih you can edit on your system.

after that run::

  atlas run all




.. 2. Download all databases first
.. -------------------------------
..
.. May be you want to make sure that all databases are downloaded correctly. Simply run::
..
..     atlas download --db-dir path/to/databases
..
.. To reassure you, most of the databases are md5 checked. The downloads use approximately 30 GB of disk space.

.. 3. Test installation
.. --------------------
..
.. Use our example_data on the GitHub repo. The first time you run atlas, it installs all dependencies.
.. It needs therefore an internet connection and some time.

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
                                    [default: /Users/silas/Documents/GitHub/atla
                                    s/databases]
    -w, --working-dir PATH          location to run atlas
    --assembler [megahit|spades]    assembler  [default: spades]
    --data-type [metagenome|metatranscriptome]
                                    sample data type  [default: metagenome]
    --interleaved-fastq             fastq files are paired-end in one files
                                    (interleaved)
    --threads INTEGER               number of threads to use per multi-threaded
                                    job
    --skip-qc                       Skip QC, if reads are already pre-processed
    -h, --help                      Show this message and exit.


This command creates a ``samples.tsv`` and a ``config.yaml`` in the working directory.

Have a look at them with a normal text editor and check if the samples names are inferred correctly.
Samples should be alphanumeric names and cam be dash delimited. Underscores should be fine too.
See the  :download:`example sample table <../reports/samples.tsv>`



The ``BinGroup`` parameter is used during the genomic binning.
In short: all samples in which you expect the same strain to
be found should belong to the same group,
e.g. all metagenome samples from mice in the same cage.
If you want to use :ref:`long reads <longreads>` for a hybrid assembly, you can also specify them in the sample table.


You should also check the ``config.yaml`` file, especially:


- You may want to add ad :ref:`host genomes <contaminants>` to be removed.
- You may want to change the resources configuration, depending on the system you run atlas on.


Details about the parameters can be found in the section :ref:`Configuration`

atlas run
----------

::

  Usage: atlas run [OPTIONS]
                   [[qc|assembly|binning|genomes|genecatalog|None|all]]
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
    --profile TEXT          snakemake profile e.g. for cluster execution.
    -n, --dryrun            Test execution.  [default: False]
    -h, --help              Show this message and exit.


``atlas run`` need to know the working directory with a ``samples.tsv`` inside it.

Take note of the ``--dryrun`` parameter, see the section :ref:`snakemake` for other handy snakemake arguments.

We recommend to use atlas on a :ref:`cluster` system, which can be set up in a view more commands.
