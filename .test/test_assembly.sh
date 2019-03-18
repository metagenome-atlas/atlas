#! /bin/bash
set -euo pipefail





atlas --version




databaseDir=".test/databases"
WD='.test/Test_assembly'
reads_dir=".test/reads/stub"

ressource_args=" --config java_mem=4 assembly_mem=4"


# gen randomreads
#very low number only for assembly
snakemake -s atlas/rules/testing.smk -d $reads_dir --config reads=1000


rm -f $WD/samples.tsv
#
atlas init --db-dir $databaseDir --threads 3  -w $WD $reads_dir


atlas run -w $WD qc $ressource_args $@

atlas run assembly -w $WD $ressource_args assembler=spades $@

rm -rf $WD/*/assembly
atlas run assembly -w $WD $ressource_args assembler=megahit $@
