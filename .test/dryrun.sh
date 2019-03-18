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
snakemake -s atlas/rules/testing.smk -d $reads_dir --config reads=1000



atlas init --db-dir $databaseDir --threads 3 -w $WD $reads_dir


atlas download --db-dir $databaseDir -n

for w in qc assembly genomes genecatalog ; do

  echo "
        Dryrun Workflow $w
      "

  atlas run $w -w $WD --dryrun $@

done
#
atlas run -w $WD --dryrun $@
#
