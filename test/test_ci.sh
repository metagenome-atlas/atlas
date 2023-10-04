#! /bin/bash

set -exuo pipefail





atlas --version

get test reads
wget https://zenodo.org/record/3992790/files/test_reads.tar.gz
tar -xzf test_reads.tar.gz 


ls -l test_reads

databaseDir="databases"
WD='test_ci'
reads_dir="test_reads" #"example_data/reads/test"


rm -f $WD/samples.tsv
#
atlas init $reads_dir --db-dir  $databaseDir  -w $WD #--interleaved-fastq

atlas run None screen -w $WD qc  $@

echo "\n\nFinished screen\n\n"

atlas run -w $WD qc  $@

echo "\n\nFinished qc\n\n"


atlas run assembly -w $WD $@

echo "\n\nFinished assembly\n\n"

atlas run binning -w $WD $@

echo "\n\nFinished binning\n\n"

atlas run genecatalog --omit-from combine_egg_nogg_annotations combine_dram_genecatalog_annotations -w $WD $@

echo "\n\nFinished genecatalog\n\n"

# atlas run genomes -w $WD $@

# echo "\n\nFinished genomes\n\n"

# atlas run all -w $WD $@

# echo "\n\nFinished all\n\n"



