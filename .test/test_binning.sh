#! /bin/bash
set -euo pipefail

atlas --version

databaseDir=".test/databases"
WD='example_data/binning'
reads_dir="example_data/reads/test"


rm -f $WD/samples.tsv $WD/finished_assembly
touch $WD/finished_assembly

#
atlas init --db-dir $databaseDir --threads 4  -w $WD --skip-qc $reads_dir

atlas run None "reports/bin_report_DASTool.html" --config interleaved_fastqs=True -w $WD $@

# genomes need databases
genome_files="genomes/clustering/allbins2genome.tsv genomes/counts/raw_counts_genomes.tsv genomes/tree/checkm.nwk"
atlas run None $genome_files --config interleaved_fastqs=True -w $WD $@
