#! /bin/bash
set -euo pipefail





atlas --version




databaseDir=".test/databases"
WD='.test/WD'


# gen randomreads
#very low number only for assembly
snakemake -s atlas/rules/testing.smk -d .test/reads --config reads=5000


rm -fr $WD
#
atlas init --db-dir $databaseDir --threads 3  -w $WD .test/reads


atlas run -w $WD qc $@

atlas run assembly -w $WD $@

#rm -rf $WD/*/assembly
#atlas run assembly -w $WD --config assembler=megahit $@
