.. _conda: http://anaconda.org/
.. _mamba: https://github.com/TheSnakePit/mamba

Getting Started
***************

Setup
=====

Conda package manager
---------------------

Atlas has **one dependency**: conda_. All databases and other dependencies are installed **on the fly**.
Atlas is based on snakemake which allows to run steps of the workflow in parallel on a cluster.

If you want to try atlas and have a linux computer (OSX may also work), you can use our `example data`_ for testing.

For real metagenomic data atlas should be run on a _linux_ sytem, with enough memory (min ~50GB but assembly usually requires 250GB).



You need to install `anaconda <http://anaconda.org/>`_ or miniconda. If you haven't done it already you need to configure conda with the bioconda-channel and the conda-forge channel. This are sources for packages beyond the default one.::

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

The order is important by the way.

Install mamba
-------------

Conda can be a bit slow because there are so many packages. A good way around this is to use mamba_ (another snake).::

    conda install mamba


From now on you can replace ``conda install`` with ``mamba install`` and see how much faster this snake is.

Install metagenome-atlas
------------------------

We recommend you to install metagenome-atlas into a conda environment e.g. named ``atlasenv``::

    mamba create -y -n atlasenv metagenome-atlas
    source activate atlasenv



Install metagenome-atlas from GitHub
------------------------------------

Alternatively you can install metagenome Atlas directly form GitHub. This allows you to access versions that are not yet in the conda release, e.g. versions that are still in development.
::


    git clone https://github.com/metagenome-atlas/atlas.git
    cd atlas
    
    # optional change to different branch
    # git checkout branchname

    # create dependencies for atlas
    mamba env create -n atlas-dev --file atlasenv.yml
    conda activate atlas-dev

    # install atlas version. Changes in this files are directly available in the atlas dev version
    pip install --editable .
    cd ..





.. _`example data`:

Example Data
============

If you want to test atlas on a small example data here is a two sample, three genome minimal metagenome dataset,
to test atlas. Even when atlas will run faster on the test data,
it will anyway download all the databases and requirements, for the a complete run,
which can take a certain amount of time and especially disk space (>100Gb).

The database dir of the test run should be the same as for the later atlas executions.

The example data can be downloaded as following::

  wget https://zenodo.org/record/3992790/files/test_reads.tar.gz
  tar -xzf test_reads.tar.gz



Usage
=====

Start a new project
-------------------

Let's apply atlas on your data or on our `example data`_::

  atlas init --db-dir databases path/to/fastq

This command creates a ``samples.tsv`` and a ``config.yaml`` in the working directory.

Have a look at them with a normal text editor and check if the samples names are inferred correctly.
Samples should be alphanumeric names and cam be dash delimited. Underscores should be fine too.
See the  :download:`example sample table <../reports/samples.tsv>`

The ``BinGroup`` parameter is used during the genomic binning.
In short: all samples in which you expect the same strain to
be found should belong to the same group,
e.g. all metagenome samples from mice in the same cage or location.
If you want to use :ref:`long reads <longreads>` for a hybrid assembly, you can also specify them in the sample table.


You should also check the ``config.yaml`` file, especially:


- You may want to add ad :ref:`host genomes <contaminants>` to be removed.
- You may want to change the resources configuration, depending on the system you run atlas on.
Details about the parameters can be found in the section :ref:`Configuration`

Keep in mind that all databases are installed in the directory specified with ``--db-dir`` so choose it wisely.


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



Run atlas
---------

::

  atlas run all


``atlas run`` need to know the working directory with a ``samples.tsv`` inside it.

Take note of the ``--dryrun`` parameter, see the section :ref:`snakemake` for other handy snakemake arguments.

We recommend to use atlas on a :ref:`cluster` system, which can be set up in a view more commands.


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



Execue Atlas
************

.. _`snakemake profile`: https://github.com/metagenome-atlas/clusterprofile

.. _cluster:

Cluster execution
=================

Automatic submitting to cluster systems
---------------------------------------

Thanks to the underlying snakemake Atlas can submit parts of the pipeline automatically to a cluster system and define the appropriate resources. If one job has finished it launches the next one.
This allows you use the full capacity of your cluster system. You even need to pay attention not to spam the other users of the cluster.




Thanks to the underlying snakemake system, atlas can submit parts of the pipeline  to clusters and cloud systems. Instead of running all steps of the pipeline in one cluster job, atlas can automatically submit each step to your cluster system, specifying the necessary threads, memory, and runtime, based on the values in the config file. Atlas periodically checks the status of each cluster job and can re-run failed jobs or continue with other jobs.

See atlas scheduling jobs on a cluster in action `<https://asciinema.org/a/337467>`_.

If you have a common cluster system (Slurm, LSF, PBS ...) we have an easy set up (see below). Otherwise, if you have a different cluster system, file a GitHub issue (feature request) so we can help you bring the magic of atlas to your cluster system.
For more information about cluster- and cloud submission, have a look at the `snakemake cluster docs <https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html>`_.

Set up of cluster execution
---------------------------

You need cookiecutter to be installed, which comes with atlas

Then run::

    cookiecutter --output-dir ~/.config/snakemake https://github.com/metagenome-atlas/clusterprofile.git

This opens a interactive shell dialog and ask you for the name of the profile and your cluster system.
We recommend you keep the default name ``cluster``. The profile was tested on ``slurm``, ``lsf`` and ``pbs``.

The resources (threads, memory and time) are defined in the atlas config file (hours and GB).

If you need to specify **queues or accounts** you can do this for all rules or for specific rules in the ``~/.config/snakemake/cluster/cluster_config.yaml``. In addition, using this file you can overwrite the resources defined  in the config file.

Example for ``cluster_config.yaml`` with queues defined::


  __default__:
  # default parameter for all rules
    queue: normal
    nodes: 1


  # The following rules in atlas need need more time/memory.
  # If you need to submit them to different queues you can configure this as outlined.

  run_megahit:
    queue: bigmem
  run_spades:
    queue: bigmem

  This rules can take longer
  run_checkm_lineage_wf:
    queue: long



Now, you can run atlas on a cluster with::

    atlas run <options> --profile cluster


As the whole pipeline can take several days, I usually run this command in a screen on the head node, even when system administrators don't normally like that. On the head node atlas only schedules the jobs and combines tables, so it doesn't use many resources. You can also submit the atlas command as a long lasting job.

 .. The mapping between  resources and cluster are defined in the ``~/.config/snakemake/cluster/key_mapping.yaml``.




If a job fails, you will find the "external jobid" in the error message.
You can investigate the job via this ID.


The atlas argument ``--jobs`` now becomes the number of jobs simultaneously submitted to the cluster system. You can set this as high as 99 if your colleagues don't mind you over-using the cluster system.


.. _local:
Single machine execution
========================

If you cannot use the  :ref:`automatic scheduling <cluster>` you can still try to use atlas on a single machine (local execution) with a lot of memory and threads ideally. In this case I recommend you the following options. The same applies if you submit a single job to a cluster running atlas.

In theory you don't need to adapt the parameters in the config file. However you should tell atlas how many threads and how much memory (GB) you have available on our system so Atlas can take this into account.

For local execution the ``--jobs`` command line arguments defines the number of threads used in total. Set it to the number of processors available on your machine.  If you have less core available than specified in the config file. The jobs are downscaled. If you have more Atlas tries to start multiple jobs, to optimally use the cores on you machine. The same applies for the memory.

For example on a machine with 16 processors and 250GB memory you might want to run::

  atlas run all --resources mem=245 --jobs 16

The whole pipeline can take more than a day. If for any reason the pipeline stops you can just rerun the same command after having inspected the error.


Cloud execution
===============

Atlas, like any other snakemake pipeline can  also easily be submitted to cloud systems. I suggest looking at the `snakemake doc <https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html>`_. Keep in mind any snakemake comand line argument can just be appended to the atlas command.



.. _snakemake:

Useful command line options
===========================

Atlas builds on snakemake. We designed the command line interface in a way that additional snakemake arguments can be added to an atlas run call.

For instance the ``--profile`` used for cluster execution. Other handy snakemake command line arguments include.

 ``--keep-going``, which  allows atlas in the case of a failed job to continue with independent steps.

For a full list of snakemake arguments see the `snakemake doc <https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options>`_.
