#! /bin/bash
set -euo pipefail





atlas --version




databaseDir=".test/databases"
WD='.test/WD'

rm -fr $WD
#
atlas init --db-dir $databaseDir --threads 3  -w $WD example_data


atlas run -w $WD qc $@

atlas run assembly -w $WD $@

atlas run assembly -w $WD $@

echo "copy qc reads and assemble"

rm -fr WD2
mkdir -p WD2/qcreads
cp WD/*/sequence_quality_control/*_QC_R?.fastq.gz WD2/qcreads


atlas init --db-dir $databaseDir --threads 3 --assembler megahit --skip-qc -w WD2 WD2/qcreads

atlas run -w WD2 assembly $snakemake_args $@
