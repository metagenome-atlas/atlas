.. _`snakemake profile`: https://github.com/metagenome-atlas/clusterprofile

.. _cluster:

Cluster execution
=================

Thank the underlying, Snakemake system, atlas can be executed on virtually all clusters and cloud systems. Instead of running all steps of the pipeline in one cluster job atlas can automatically submit each step to your cluster system, specifying the necessary threads, memory, and runtime, based on the values in the config file. Atlas periodically checks the status of each cluster job and can re-run failed jobs or continue with other jobs.


If you have a common cluster system (Slurm, LSF, PBS ...) we have an easy set up (see below). Otherwise, if you have an other cluster system, fill in an GitHub issue (feature request) so we can help you bringing the magic of atlas to your cluster system.
For more information about cluster- and cloud submission, have a look at the `snakemake doc <https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html>`_.

Set up of cluster execution
---------------------------

You need cookiecutter to be installed, which comes with atlas

Then run::

    cookiecutter --output-dir ~/.config/snakemake https://github.com/metagenome-atlas/clusterprofile.git

This opens a interactive shell dialog and ask you for the name of the profile and your cluster system.
We recommend you to keep the default name ``cluster``. The profile supports ``slurm``,``lsf`` and ``pbs``.

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


 .. The mapping between  resources and cluster are defined in the ``~/.config/snakemake/cluster/key_mapping.yaml``.




If a job fails, you will find the "external jobid" in the error message.
You can investigate the job via this ID.


Useful command line options
----------------------------

The atlas argument ``--jobs`` now becomes the number of jobs simultaneously submitted to the cluster system. You can set this as high as 99 if you don't have a problem with your colleges of over-using the cluster system.

In the case of an failed job, ``--keep-going`` (default false)  allows atlas to continue with independent steps.


Cloud execution
===============

Atlas, like any other snakemake pipeline can thanks also easily be submitted to cloud systems. I let look at the `snakemake doc <https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html>`_. Keep in mind any snakemake comand line argument can just be appended to the atlas command.

.. _local:

Local execution
===============
The number of threads used **for each step** can be configured in the config file::

  threads: 8
  assembly_threads: 8

For local execution the ``--jobs`` command line arguments defines the number of threads used in total and by default is set to the number of processor of your system. If you have less core available than specified in the config file. The jobs are downscaled. If you have more Atlas tries to start multiple jobs, to optimally use the cores on you machine.
If you don't want that Atlas uses all the available cores then you can specify this with the ``--jobs`` command line argument.
