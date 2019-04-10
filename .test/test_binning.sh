#! /bin/bash
set -euo pipefail





atlas --version




databaseDir=".test/databases"
WD='.test/Test_binning'
reads_dir=".test/reads/binning"

ressource_args=" --config java_mem=10 assembly_mem=10"


# gen randomreads
snakemake -s atlas/rules/testing.smk -d $reads_dir --config reads=50000


rm -f $WD/samples.tsv
#
atlas init --db-dir $databaseDir --threads 3  -w $WD $reads_dir

atlas run binning -w $WD $ressource_args assembler=spades final_binner=metabat $@

# genomes need databases
atlas run genomes -w $WD $ressource_args assembler=spades final_binner=metabat --omit-from get_genome_for_cat $@

#atlas run all -w $WD $ressource_args assembler=spades final_binner=metabat $@
