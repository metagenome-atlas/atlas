Defining Samples
================

Samples are defined with a name, file path(s) and the type of data. A single
file path is interpreted as interleaved paired-end reads or single-end, while
two paths must include full paths to R1 and R2.

Sample names must be unique and not contain spaces or underscores (dashes are accepted).

For 'type', the value can be either 'metagenome' or 'metatranscriptome'. The
default if not specified is 'metagenome'.


Interleaved PE input::

    samples:
        sample-1:
            fastq: /data/sample-1_pe.fastq.gz
            type: metagenome


Paired-end as separate files::

    samples:
        sample-1:
            fastq:
                - /data/sample-1_R1.fastq.gz
                - /data/sample-1_R2.fastq.gz
            type: metagenome

The '-' is required if multiple fastq file paths need to be specified.

Data is assumed to be paired-end unless stated otherwise. If your data is
single-end sequence data, specify 'paired' as ``false``::

    samples:
        sample-1:
            fastq: /data/sample-1_pe.fastq.gz
            type: metagenome
            paired: false
