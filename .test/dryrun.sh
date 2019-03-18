#! /bin/bash
set -euo pipefail


rm -rf WD



atlas --version
atlas run --help



databaseDir=".test/databases"
WD='.test/Dryrun'
reads_dir=".test/reads/stub"

rm -fr $WD
#
# gen randomreads
#very low number only for assembly
#snakemake -s atlas/rules/testing.smk -d $reads_dir --config reads=1000



atlas init --db-dir $databaseDir --threads 3 -w $WD $reads_dir


atlas download --db-dir $databaseDir -n


atlas init --db-dir $databaseDir --threads 3 -w $WD --assembler spades example_data


for w in qc assembly genomes genecatalog ; do

  echo "
        Dryrun Workflow $w
      "

  atlas run $w -w $WD --dryrun $@

done
#
atlas run -w $WD --dryrun $@
#

# skip QC

WD=${WD}_noQC
rm -fr $WD
atlas init --db-dir $databaseDir --threads 3 --skip-qc -w $WD --assembler megahit example_data


for w in assembly genomes genecatalog ; do

  echo "
        Dryrun Workflow $w
      "

  atlas run $w -w $WD --dryrun $@

done
#
atlas run -w $WD --dryrun $@
