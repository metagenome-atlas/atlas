#! /bin/bash
set -euo pipefail


rm -rf WD



atlas --version
atlas run --help


databaseDir=".test/databases"
WD='.test/WD'

rm -fr $WD
#


atlas download --db-dir $databaseDir -n


atlas init --db-dir $databaseDir --threads 3 -w $WD example_data


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

atlas init --db-dir $databaseDir --threads 3 --skip-qc -w $WD example_data


for w in assembly genomes genecatalog ; do

  echo "
        Dryrun Workflow $w
      "

  atlas run $w -w $WD --dryrun $@

done
#
atlas run -w $WD --dryrun $@
