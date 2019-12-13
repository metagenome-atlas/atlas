#! /bin/bash
set -euo pipefail

atlas --version

databaseDir="databases"
WD='example_data/binning'
reads_dir="example_data/reads/test"
config="--config interleaved_fastqs=True threads=2 mem=4 java_mem=4 "

rm -f $WD/samples.tsv $WD/finished_assembly
touch $WD/finished_assembly

#
atlas init --db-dir $databaseDir --threads 4  -w $WD --skip-qc $reads_dir

atlas run None binning $config -w $WD $@ --omit-from get_bins

echo "Copy checkm files\n\n\n"
pwd
for s in "sample1 sample2" ;
do
  dest_dir=$WD/$s/binning/DASTool/checkm/
  mkdir -p $dest_dir
  tree $WD/checkm_results/
  mv $WD/checkm_results/$s/* $dest_dir
done

atlas run None binning $config -w $WD $@

# genomes need databases
atlas run None genomes/clustering/allbins2genome.tsv  $config -w $WD $@
