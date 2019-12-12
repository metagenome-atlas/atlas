#! /bin/bash
set -euo pipefail

atlas --version

databaseDir="databases"
WD='example_data/binning'
reads_dir="example_data/reads/test"
config="--config interleaved_fastqs=True final_binner=maxbin threads=2 mem=4 java_mem=4"

rm -f $WD/samples.tsv $WD/finished_assembly
touch $WD/finished_assembly

#
atlas init --db-dir $databaseDir --threads 4  -w $WD --skip-qc $reads_dir

atlas run None "reports/bin_report_metabat.html" $config -w $WD $@

# genomes need databases
genome_files="genomes/clustering/allbins2genome.tsv genomes/counts/raw_counts_genomes.tsv genomes/tree/checkm.nwk"
atlas run None $genome_files $config -w $WD $@
