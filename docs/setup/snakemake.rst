Snakemake
=========

Command-line execution of Snakemake is handled by ATLAS, however, a user
can still opt to use Snakemake outside of the ATLAS command line interface.

When executing a command from ``atlas``, the first line of output will include
the snakemake command. For example, if we are getting started and need to
download the databases, we execute::

    atlas download --jobs 1 --out-dir metaomics/references

Here is the abbreviated output from the above command::

    [2017-06-30 11:48 INFO] Executing: snakemake -s /people/brow015/anaconda3/lib/python3.5/site-packages/atlas/Snakefile -d /pic/projects/mint/metaomics -p -j 1 --nolock --rerun-incomplete --config db_dir='/pic/projects/mint/metaomics/references' workflow=download --
    Provided cores: 1
    Rules claiming more threads will be scaled down.
    Job counts:
    	count	jobs
    	1	all
    	6	transfer_files
    	7

    rule transfer_files:
        output: /pic/projects/mint/metaomics/references/silva_rfam_all_rRNAs.fa
        jobid: 1
        wildcards: filename=silva_rfam_all_rRNAs.fa

    curl 'https://zenodo.org/record/804435/files/silva_rfam_all_rRNAs.fa' -s > /pic/projects/mint/metaomics/references/silva_rfam_all_rRNAs.fa

    ...
    6 of 7 steps (86%) done

    localrule all:
        input: /pic/projects/mint/metaomics/references/adapters.fa, /pic/projects/mint/metaomics/references/phiX174_virus.fa, /pic/projects/mint/metaomics/references/silva_rfam_all_rRNAs.fa, /pic/projects/mint/metaomics/references/refseq.tree, /pic/projects/mint/metaomics/references/refseq.dmnd, /pic/projects/mint/metaomics/references/refseq.db
        jobid: 0

    Finished job 0.
    7 of 7 steps (100%) done
    All databases have downloaded and validated successfully.
    When generating your configuration file, use '--database-dir /pic/projects/mint/metaomics/references'

We can see that the workflow being utilized is::

    /people/brow015/anaconda3/lib/python3.5/site-packages/atlas/Snakefile

So any later or more complex use cases can simply call ``snakemake`` and
specify ``--snakefile`` for their local instance.


Provenance
----------

Extra arguments on the command line are passed directly into the snakemake
call, so even within ``atlas assemble`` we can do things like::

    atlas assemble --jobs 24 --out-dir test-dir config.yaml --summary

This results in the call of::

    snakemake -s /anaconda3/lib/python3.5/site-packages/atlas/Snakefile \
        -d /pic/projects/mint/metaomics/results/subset \
        -p -j 24 --rerun-incomplete \
        --configfile config.yaml \
        --nolock \
        --use-conda \
        --config workflow=complete \
        --summary

The output gives details per output file, which rule created the file, the
creation date, and other information that is relevant to the file's creation.

Snakemake's ``--detailed-summary`` adds columns for input file as well as the
shell command that was used.
