.. _execution_system:

Execution of Atlas
==================
otherwise use multiple cores. The number of threads used **for each step** can be configured in the config file::

  threads: 8
  assembly_threads: 8

When you execute Atlas it checks how many cores are available. If you have less core available than specified in the config file.
The jobs are downscaled. If you have more Atlas tries to start multiple jobs, to optimally use the cores on you machine.
If you don't want that Atlas uses all the available cores then you can specify this with the ``--jobs`` command line argument.


Cluster execution
-----------------

We recommend you to execute Atlas on a cluster system. Thank to the underlying, Snakemake workflow, Atlas can be executed on virtually all cluster systems.
Each job get submitted to your cluster system with the memory (in GB) and threads specified in the config file. The ``--jobs`` command line argument now defines
*how many jobs you want to run simultaneously on the cluster system.* The default is still the number of cores on your computer but you can set it higher ``--jobs 99``.

Have a look at the snakemake doc how to execute snakemake pipelines on clusters. If you have a common cluster system (Slurm, LSF, toque ...) you can use our `snakemake profile`_.
Note, in atlas memory is specified in GB.

.. _`snakemake profile`: https://github.com/metagenome-atlas/generic
