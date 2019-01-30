The Atlas pipeline
==================

|scheme|

.. |scheme| image:: ../resources/images/ATLAS_scheme.png
  :alt: Atlas is a workflow for assembly and binning of metagenomic reads


Quality control
---------------

::
  atlas run qc


Runs quality control of single or paired end reads.

Main output files:

  - `reports/QC_report.html`_
  - Different statis in ``stats``


.. _reports/QC_report.html: reports/QC_report.html

Per sample:

  - {sample}/sequence_quality_control/{sample}_QC_{fraction}.fastq.gz
  - Various quality stats in {sample}/sequence_quality_control/read_stats

Fractions:
``````````
When the input was paired end, we will put out three the reads in three fractions R1,R2 and se
The se are the paired end reads which lost their mate during the filtering.
The se are seamlessly integrated in the next steps.


Assembly
---------------

Genomes
---------------


Gene Catalog
---------------
