#!/bin/bash
set -euo pipefail


NThreads=2
MaxMem=3

atlas --version
atlas run --help


databaseDir="test/databases"
WD='test/Dryrun'
reads_dir='test/reads/empty'
snakemake_args=" --quiet rules $@ --dryrun " 
test_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"



create_reads_dir() {

    local reads_dir="$1"
    local N=$2

echo "touch reads dir: $reads_dir"

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

# need at least 10 samples for cobinning

create_reads_dir $reads_dir 10 





rm -fr $WD

echo "Atlas download"
atlas download --db-dir $databaseDir -n

echo "Init"
atlas init --db-dir $databaseDir --threads=$NThreads -w $WD $reads_dir




echo "Dryrun all"
atlas run all -w $WD  $snakemake_args

echo "Dryrun strains"
atlas run genomes strains -w $WD $snakemake_args


for binner in metabat SemiBin vamb DASTool ; do

  echo "
        Dryrun Binner $binner
      "

  atlas run binning -w $WD --config final_binner=$binner $snakemake_args

done


#

echo "
      Dryrun with skip QC and megahit
    "
#

rm -fr $WD

WD=${WD}/noQC
rm -fr $WD

atlas init --db-dir $databaseDir --skip-qc -w $WD --assembler megahit $reads_dir

atlas run all -w $WD $snakemake_args


echo "
      execution with profile
    "

  mkdir -p $WD/local
  printf 'cores: 2\n' > $WD/local/config.yaml

  atlas run qc -w $WD  --profile $WD/local $snakemake_args


# clean up
rm -rf $WD $reads_dir






echo " 
      test with external genomes
      "

bash $test_script_dir/test_external_genomes.sh $snakemake_args



echo " 
      test init with different samples
      "

bash $test_script_dir/test_init_many_samples.sh $snakemake_args