#!/bin/bash
set -euo pipefail


NThreads=2
MaxMem=3

atlas --version
atlas run --help


databaseDir="test/databases"
WD='test/Dryrun'
reads_dir='test/reads/empty'


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

echo "Atlas download"
atlas download --db-dir $databaseDir -n

echo "Init"
atlas init --db-dir $databaseDir --threads=$NThreads -w $WD $reads_dir




echo "Dryrun all"
atlas run all -w $WD --max-mem $MaxMem --jobs $NThreads --dryrun $@

echo "Dryrun strains"
atlas run genomes strains -w $WD --max-mem $MaxMem --jobs $NThreads --dryrun $@


for binner in SemiBin vamb DASTool ; do

  echo "
        Dryrun Binner $binner
      "

  atlas run genomes -w $WD --config final_binner=$binner --dryrun $@

done

atlas run quantify_genomes -w $WD  --dryrun $@
#

echo "
      Dryrun with skip QC and megahit
    "
#
WD=${WD}/noQC
rm -fr $WD
atlas init --db-dir $databaseDir --threads=$NThreads --skip-qc -w $WD --assembler megahit $reads_dir

atlas run all -w $WD --dryrun $@


echo

echo "
      execution with profile
    "

  mkdir -p $WD/local
  printf 'cores: 1\n' > $WD/local/config.yaml

  atlas run qc -w $WD --dryrun --max-mem $MaxMem --jobs $NThreads --profile $WD/local $@ 
