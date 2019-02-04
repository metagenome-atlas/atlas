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
  - Different stats in ``stats``


.. _reports/QC_report.html: reports/QC_report.html

.. raw:: html
    :file: ../reports/QC_report.html

Per sample:

  - {sample}/sequence_quality_control/{sample}_QC_{fraction}.fastq.gz
  - Various quality stats in {sample}/sequence_quality_control/read_stats

.. _fractions:

Fractions:
``````````
When the input was paired end, we will put out three the reads in three fractions R1,R2 and se
The se are the paired end reads which lost their mate during the filtering.
The se are seamlessly integrated in the next steps.


Assembly
---------------

Besides the `reports/assembly_report.html`_ this rule outputs the following files per sample:

  -  "{sample}/{sample}_contigs.fasta",
  - "{sample}/sequence_alignment/{sample}.bam",
  -  "{sample}/assembly/contig_stats/postfilter_coverage_stats.txt",
  -      "{sample}/assembly/contig_stats/prefilter_contig_stats.txt",
  -      "{sample}/assembly/contig_stats/final_contig_stats.txt"


.. _reports/assembly_report.html:



Genomes
---------------


Gene Catalog
---------------
