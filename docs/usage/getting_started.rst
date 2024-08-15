.. _conda: http://anaconda.org/
.. _mamba: https://github.com/TheSnakePit/mamba

Getting Started
***************

Setup
=====

Conda package manager
---------------------

Atlas has **one dependency**: conda_. All databases and other dependencies are installed **on the fly**.
Atlas is based on snakemake, which allows to run steps of the workflow in parallel on a cluster.

If you want to try atlas and have a linux computer (OSX may also work), you can use our `example data`_ for testing.

For real metagenomic data atlas should be run on a _linux_ system, with enough memory (min ~50GB but assembly usually requires 250GB).



You need to install `anaconda <http://anaconda.org/>`_ or miniconda. 
If you haven't done it already, you need to configure conda with the bioconda-channel and the conda-forge channel. This are sources for packages beyond the default one.
Setting strict channel priority can prevent quite some annoyances.

.. code-block:: bash
    conda config --set channel_priority strict
    conda config --add channels bioconda
    conda config --add channels conda-forge

The order is important by the way.

Install mamba
-------------

Conda can be a bit slow because there are so many packages. A good way around this is to use mamba_ (another snake).::

    conda install mamba


From now on, you can replace ``conda install`` with ``mamba install`` and see how much faster this snake is.

Install metagenome-atlas
------------------------

We recommend to install metagenome-atlas into a conda environment e.g. named ``atlasenv``. 
We also recommend to specify the latest version of metagenome-atlas.  

.. code-block:: bash

    mamba create -y -n atlasenv metagenome-atlas={latest_version}
    source activate atlasenv

where `{latest_version}` should be replaced by 

.. image:: https://anaconda.org/bioconda/metagenome-atlas/badges/version.svg
    :target: https://anaconda.org/bioconda/metagenome-atlas




Install metagenome-atlas from GitHub
------------------------------------

Alternatively, you can install metagenome Atlas directly from GitHub. This allows you to access versions that are not yet in the conda release, e.g. versions that are still in development.

.. code-block:: bash

    git clone https://github.com/metagenome-atlas/atlas.git
    cd atlas

    # optional change to different branch
    # git checkout branchname

    # create dependencies for atlas
    mamba env create -n atlas-dev --file atlasenv.yml
    conda activate atlas-dev

    # install atlas version. Changes in the files are directly available in the atlas dev version
    pip install --editable .
    cd ..





.. _`example data`:

Example Data
============

If you want to test atlas on a small example data, here is a two sample, three genome minimal metagenome dataset,
to test atlas. Even when atlas will run faster on the test data,
it will anyway download all the databases and requirements, for a complete run,
which can take a certain amount of time and especially disk space (>100Gb).

The database dir of the test run should be the same as for the later atlas executions.

The example data can be downloaded as following

.. code-block:: bash

  wget https://zenodo.org/record/3992790/files/test_reads.tar.gz
  tar -xzf test_reads.tar.gz



Usage
=====

Start a new project
-------------------

Let's apply atlas on your data or on our `example data`_::

  atlas init --db-dir databases path/to/fastq_files

This command parses the folder for fastq files (extension ``.fastq(.gz)`` or ``.fq(.gz)`` , gzipped or not). fastq files can be arranged in subfolders, in which case the subfolder name will be used as a sample name. If you have paired-end reads the files are usually distinguishable by ``_R1/_R2`` or simple ``_1/_2`` in the file names. Atlas searches for these patterns and lists the paired-end files for each sample.

The command creates a ``samples.tsv`` and a ``config.yaml`` in the working directory.

Have a look at them with a normal text editor and check if the sample names are inferred correctly. The sample names are used for the naming of contigs, genes, and genomes. Therefore, the sample names should consist only of digits and letters and start with a letter (Even though one ``-`` is allowed). Atlas tries to simplify the file name to obtain unique sample names, if it doesn't succeed it simply puts S1, S2, ... as sample names.


See the  :download:`example sample table <../reports/samples.tsv>`

The ``BinGroup`` parameter is used during the genomic binning.
In short: If you have between 5 and 150 samples the default (putting everything in one group) is fine.
If you have less than 5 samples, put every sample in an individual BinGroup and use `metabat` as final binner.
If you have more samples see the :ref:`cobinning` section for more details.

.. note:: If you want to use :ref:`long reads <longreads>` for a hybrid assembly, you can also specify them in the sample table.


You should also check the ``config.yaml`` file, especially:


- You may want to add ad :ref:`host genomes <contaminants>` to be removed.
- You may want to change the resources configuration, depending on the system you run atlas on.

Details about the parameters can be found in the section :ref:`Configuration`

Keep in mind that all databases are installed in the directory specified with ``--db-dir`` so choose it wisely.


.. code-block:: text

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



Start a new project with public data
------------------------------------

Since v2.9 atlas has possibility to start a new project from public data stored in the short read archive (SRA).

You can run ``atlas init-public <SRA_IDs>`` and specify any ids, like bioprojects, or other SRA ids. 

Atlas does the following steps:

  1. Search SRA for the corresponding sequences (Runs) and save them in the file ``SRA/RunInfo_original.tsv``. For example, if you specify a Bioproject, it fetches the information for all runs of this project. 
  2. Atlas filters the runs to contain only valid metagenome sequences. E.g. exclude singleton reads, 16S. The output will be saved in ``RunInfo.tsv``
  3. Sometimes the same Sample is sequenced on different lanes, which will result into multiple runs from the same sample. Atlas will **merge** runs from the same biosample.
  4. Prepare a sample table and a config.yaml similar to the ``atlas init`` command.


If you are not happy with the filtering atlas performs, you can go back to the ``SRA/RunInfo_original.tsv`` and create a new ``RunInfo.tsv``. 
If you then rerun ``atlas init-public continue`` it will continue from your modified RunInfo and do step 3. & 4. above. 


Limitations: For now atlas, cannot handle a mixture of paired and single end reads, so we focus primarily on the paired end. 
If you have longreads for your project, you would need to specify them yourself in the sample.tsv.

During the run, the reads are downloaded from SRA in the likely most efficient way using prefetch and parallel, fastq.gz generation. 
The download step has checkpoints, so if the pipeline gets interrupted, you can restart where you left off. 
Using the command line arguments ``--restart-times 3 and --keep-going`` You can even ask atlas to do multiple restarts before stopping. 

The downloaded reads are directly processed. However, if you only want to download the reads you can use::

  atlas run None download_sra

Example: Downloading reads from the human microbiome project2
`````````````````````````````````````````````````````````````
::

  atlas init-public --working-dir HMP2 PRJNA398089

Gives the output::
  
  [Atlas] INFO: Downloading runinfo from SRA
  [Atlas] INFO: Start with 2979 runs from 2979 samples
  [Atlas] INFO: Runs have the following values for library_source: METAGENOMIC, METATRANSCRIPTOMIC
          Select only runs library_source == METAGENOMIC, Filtered out 762 runs
  [Atlas] INFO: Runs have the following values for library_selection: PCR, RT-PCR, RANDOM
          Select only runs library_selection == RANDOM, Filtered out 879 runs
  [Atlas] INFO: Selected 1338 runs from 1338 samples
  [Atlas] INFO: Write filtered runinfo to HMP2/RunInfo.tsv
  [Atlas] INFO: Prepared sample table with 1338 samples
  [Atlas] INFO: Configuration file written to HMP2/config.yaml
          You may want to edit it using any text editor.





Run atlas
---------

::

  atlas run genomes


``atlas run`` need to know the working directory with a ``samples.tsv`` inside it.

Take note of the ``--dryrun`` parameter, see the section :ref:`snakemake` for other handy snakemake arguments.

We recommend to use atlas on a :ref:`cluster` system, which can be set up in a view more commands.


.. code-block:: text

  Usage: atlas run [OPTIONS] [qc|assembly|binning|genomes|genecatalog|None|all]
                   [SNAKEMAKE_ARGS]...

    Runs the ATLAS pipeline

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
                            submission for more details).

    --profile TEXT          snakemake profile e.g. for cluster execution.
    -n, --dryrun            Test execution.  [default: False]
    -h, --help              Show this message and exit.


Execute Atlas
************


.. _cluster:

Cluster execution
=================

Automatic submitting to cluster systems
---------------------------------------

Thanks to the underlying snakemake Atlas can submit parts of the pipeline automatically to a cluster system and define the appropriate resources. If one job has finished it launches the next one.
This allows to use the full capacity of your cluster system. You even need to pay attention not to spam the other users of the cluster.




Thanks to the underlying snakemake system, atlas can submit parts of the pipeline  to clusters and cloud systems. Instead of running all steps of the pipeline in one cluster job, atlas can automatically submit each step to your cluster system, specifying the necessary threads, memory, and runtime, based on the values in the config file. Atlas periodically checks the status of each cluster job and can re-run failed jobs or continue with other jobs.

See atlas scheduling jobs on a cluster in action `<https://asciinema.org/a/337467>`_.

If you have a common cluster system (Slurm, LSF, PBS ...) we have an easy set up (see below). Otherwise, if you have a different cluster system, file a GitHub issue (feature request) so we can help you bring the magic of atlas to your cluster system.
For more information about cluster- and cloud submission, have a look at the `snakemake cluster docs <https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html>`_.

Set up of cluster execution
---------------------------

You need cookiecutter to be installed, which comes with atlas

Then run::

    cookiecutter --output-dir ~/.config/snakemake https://github.com/metagenome-atlas/clusterprofile.git

This opens an interactive shell dialog and ask you for the name of the profile and your cluster system.
We recommend you keep the default name ``cluster``. The profile was tested on ``slurm``, ``lsf`` and ``pbs``.

The resources (threads, memory and time) are defined in the atlas config file (hours and GB).

**Specify queues and accounts**


If you have different **queues/partitions** on your cluster system you should tell atlas about them so it can *automatically choose the best queue*. Adapt the template for the queues.tsv::

  cp ~/.config/snakemake/cluster/queues.tsv.example ~/.config/snakemake/cluster/queues.tsv

Now enter the information about the queues/partitions on your particular system.


If you need to specify **accounts** or other options for one or all rules you can do this for all rules or for specific rules in the ``~/.config/snakemake/cluster/cluster_config.yaml``. In addition, using this file you can overwrite the resources defined  in the config file.

Example for ``cluster_config.yaml`` with queues defined::


  __default__:
  # default parameter for all rules
    account: project_1345
    nodes: 1



Now, you can run atlas on a cluster with::

    atlas run <options> --profile cluster


As the whole pipeline can take several days, I usually run atlas itself on a cluster in a long running queue. 

 .. The mapping between  resources and cluster are defined in the ``~/.config/snakemake/cluster/key_mapping.yaml``.


If a job fails, you will find the "external jobid" in the error message.
You can investigate the job via this ID.


The atlas argument ``--jobs`` now becomes the number of jobs simultaneously submitted to the cluster system. You can set this as high as 99 if your colleagues don't mind you over-using the cluster system.


.. _local:

Single machine execution
========================

If you don't want to use the  :ref:`automatic scheduling <cluster>` you can use atlas on a single machine (local execution) with a lot of memory and threads ideally. In this case I recommend you the following options. The same applies if you submit a single job to a cluster running atlas.

Atlas detects how many CPUs and how much memory is available on your system and it will schedule as many jobs in parallel as possible.  If you have less resources available than specified in the config file, the jobs are downscaled.

By default atlas will use all cpus and 95% of all the available memory. If you are not happy with that, or you need to specify an exact amount of memory/ cpus you can use the command line arguments ``--jobs`` and ``--max-mem`` to do so. 


Cloud execution
===============

Atlas, like any other snakemake pipeline can  also easily be submitted to cloud systems. I suggest looking at the `snakemake doc <https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html>`_. Keep in mind any snakemake command line argument can just be appended to the atlas command.



.. _snakemake:

Useful command line options
===========================

Atlas builds on snakemake. We designed the command line interface in a way that additional snakemake arguments can be added to an atlas run call.

For instance the ``--profile`` used for cluster execution. Other handy snakemake command line arguments include:

 ``--keep-going``, which allows atlas in the case of a failed job to continue with independent steps.

 ``--report``, which allows atlas to generate a user-friendly run report (e.g., by specifying ``--report report.html``). This report includes the steps used in the analysis workflow and the versions of software tools used at each step. See discussions `#523 <https://github.com/metagenome-atlas/atlas/discussions/523>`_ and `#514 <https://github.com/metagenome-atlas/atlas/discussions/514))>`_.

For a full list of snakemake arguments see the `snakemake doc <https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options>`_.
