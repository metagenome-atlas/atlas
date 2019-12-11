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
atlas run None "genomes/clustering/allbins2genome.tsv" "genomes/counts/raw_counts_genomes.tsv" --config interleaved_fastqs=True -w $WD $@

#atlas run all -w $WD $ressource_args assembler=spades final_binner=metabat $@
