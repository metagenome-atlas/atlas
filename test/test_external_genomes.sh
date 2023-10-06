#!/bin/bash
set -euo pipefail


NThreads=2
MaxMem=3





databaseDir="test/databases"
WD='test/genome_quant'
reads_dir='test/reads/empty'
snakemake_args=" --quiet rules $@ --dryrun " 

echo "touch reads dir"
mkdir -p $reads_dir
for sample in Sample1 Sample2 ; 
  do
  for fraction in R1 R2;
    do
    touch $reads_dir/${sample}_${fraction}.fastq.gz
  done
done


rm -fr $WD

# create genome dir
genome_dir=$WD/other_genomes

mkdir -p $genome_dir
for i in 1::5 ; 
do
    touch $genome_dir/Genome_$i.fasta
done

echo "Init"
atlas init --db-dir $databaseDir --skip-qc -w $WD $reads_dir

echo "Run"

atlas run quantify_genomes -w $WD --config genome_dir="other_genomes" $snakemake_args


rm -rf $WD $reads_dir