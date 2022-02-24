#!/bin/bash
set -euo pipefail


atlas --version
atlas run --help


databaseDir="test/databases"
WD='test/Dryrun'
reads_dir='test/reads/empty'


# touch reads dir
mkdir -p $reads_dir
for sample in Sample1 Sample2 ; 
  do
  for fraction in R1 R2;
    do
    touch $reads_dir/${sample}_${fraction}.fastq.gz
  done
done


rm -fr $WD


atlas download --db-dir $databaseDir -n


atlas init --db-dir $databaseDir --threads 3 -w $WD $reads_dir




echo "Dryrun all"
atlas run all -w $WD --dryrun $@

for binner in SemiBin vamb DASTool ; do

  echo "
        Dryrun Binner $binner
      "

  atlas run genomes -w $WD --config final_binner=$binner --dryrun $@

done


#

echo "Dryrun with skip QC and megahit"
#
WD=${WD}/noQC
rm -fr $WD
atlas init --db-dir $databaseDir --threads 3 --skip-qc -w $WD --assembler megahit $reads_dir

atlas run all -w $WD --dryrun $@


echo "local execution"
  atlas run qc -w $WD --dryrun --jobs 3 --max-mem 2 $@