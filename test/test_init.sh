#!/bin/bash
set -euo pipefail


NThreads=2
MaxMem=3

atlas --version
atlas run --help


databaseDir="test/databases"




create_reads_dir() {

    local reads_dir="$12"
    local N=$2

echo "touch reads dir"

rm -rf $reads_dir
mkdir -p $reads_dir

for (( i=1; i<=$N; i++ )); do
    sample="Sample$i"
    
  for fraction in R1 R2;
    do
    touch $reads_dir/${sample}_${fraction}.fastq.gz
  done
done
}



echoe "test to create projects for N samples"

for N in 5 10 50 300 ; 
do

WD="test/test_init/$N"

reads_dir="test/test_init/reads_${N}_samples/"

create_reads_dir $reads_dir $N

atlas init --db-dir $databaseDir -w $WD $reads_dir

done
