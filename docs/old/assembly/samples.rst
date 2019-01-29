Defining Samples
================

Samples are defined with a name, file path(s) and the type of data. A single
file path is interpreted as interleaved paired-end reads or single-end, while
two paths must include full paths to R1 and R2.

Sample names must be unique and not contain spaces or underscores (dashes are
accepted).

For ``type``, the value can be either 'metagenome' or 'metatranscriptome'. If
neither is specified, the default is 'metagenome'.


Interleaved Input
-----------------

A single file path is specified for ``fastq`` with ``type`` and ``paired``
also being set. In this case, ``type`` and ``paired`` are optional as we are
using the values that are equal to the defaults.

::

    samples:
        sample-1:
            fastq: /data/sample-1_pe.fastq.gz
            type: metagenome
            paired: true


Paired-end Input
----------------

In this case, we create a list using YAML_ syntax for both R1 and R2 indexes::

    samples:
        sample-1:
            fastq:
                - /data/sample-1_R1.fastq.gz
                - /data/sample-1_R2.fastq.gz
            type: metagenome
            paired: true

The '-' is required if multiple fastq file paths need to be specified.


Single-end Input
----------------

Data is assumed to be paired-end unless stated otherwise. If your data is
single-end sequence data, specify ``paired`` as ``false``::

    samples:
        sample-1:
            fastq: /data/sample-1_pe.fastq.gz
            type: metagenome
            paired: false


.. _YAML: http://www.yaml.org/
