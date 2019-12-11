#! /bin/bash
set -euo pipefail


atlas --version
atlas run --help


databaseDir=".test/databases"
WD='.test/Dryrun'
reads_dir='.test/reads/empty'

rm -fr $WD


atlas download --db-dir $databaseDir -n


atlas init --db-dir $databaseDir --threads 3 -w $WD $reads_dir


# for w in qc assembly genomes genecatalog ; do
#
#   echo "
#         Dryrun Workflow $w
#       "
#
#   atlas run $w -w $WD --dryrun $@
#
# done
#
atlas run -w $WD --dryrun $@
#

# # skip QC
#
# WD=${WD}/noQC
# rm -fr $WD
# atlas init --db-dir $databaseDir --threads 3 --skip-qc -w $WD --assembler megahit $reads_dir
#
#
# atlas run all -w $WD --dryrun $@
