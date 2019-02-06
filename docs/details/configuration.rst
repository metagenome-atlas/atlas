
.. _configuration:

Configure Atlas
===============

Example config file
-------------------


.. include:: ../../atlas/template_config.yaml
  :code:



.. _contaminants:

Remove reads from Host
----------------------

One of the most important steps in the Quality control is to remove host genome.
You can add any number of genomes to be removed.

We recommend you to use genomes where repetitive sequences are masked.
See here for more details `human genome <http://seqanswers.com/forums/archive/index.php/t-42552.html>`_.



Detailed configuration
----------------------

.. toctree::
    :maxdepth: 1

    ../advanced/qc
    ../advanced/assembly
