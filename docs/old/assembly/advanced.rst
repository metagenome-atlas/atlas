Advanced Options
================


Temporary Directory
-------------------

Some steps, like assembly, may have an optional temporary directory. If
specified in the configuration or as an environmental variable
[``TMPDIR``|``TEMP``|``TMP``], the step will use this directory rather than
the current working directory.

**Default: the working directory**

::

    temporary_directory: /scratch


Shell Command Prefix
--------------------

This allows the user to prepend something to each command. In our case, we
reserve a SLURM allocation, then submit jobs across the reservation with::

    prefix: srun --exclusive -N1 -n1 -c__threads__

And ``threads`` from the top-level of the configurations is prepended to each
multi-threaded step as the workflow is executed.

Additionally, for multi-threaded jobs, the user can set environmental variable
`SHPFXM` while for single-threaded jobs the user can set `SHPFXS` as the
command prefix.

**Default: None**
