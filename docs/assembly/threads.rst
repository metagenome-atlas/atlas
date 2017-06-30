Jobs and Threads
================

Most steps of the workflow are utilizing applications that can thread or
otherwise use multiple cores. Leaving this one below the max, in cases where
many samples are being analyzed, may be optimal as single-threaded jobs will
be processed more efficiently.

.. important::
    Threads is used per step and will most likely be a subset of ``--jobs``.
    ``--jobs`` represents the total available cores for all simultaneous steps.

**Default: 1**

::

    threads: 23

When starting your ``atlas`` command, e.g. ``atlas assemble --jobs 48 config.yaml``,
be sure to set the total thread pool to capture all available possible jobs to
be executed simultaneously. For example, if we are utilizing 3 nodes, each
with 24 cores, we would set ``threads: 24`` and execute ``atlas`` with::

    atlas assemble --jobs 72 config.yaml
