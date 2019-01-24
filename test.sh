#! /bin/bash
set -euo pipefail


#rm -rf WD



atlas --version
atlas run --help



databaseDir=".test/databases"
WD='.test/WD'

rm -f $WD/samples.tsv
#
atlas init --db-dir $databaseDir --threads 3 --assembler spades -w $WD example_data

for w in qc assembly genomes genecatalog ; do

  echo "
        Dryrun Workflow $w
      "

  atlas run all -w $WD $w --dryrun $@

done
#
atlas run -w $WD --dryrun $@
#
atlas run -w $WD qc $@



#mkdir -p WD2/qcreads
#cp WD/*/sequence_quality_control/*_QC_R?.fastq.gz WD2/qcreads

#rm -f WD2/samples.tsv

#atlas init --db-dir $databaseDir --threads 3 --assembler spades --skip-qc -w WD2 WD2/qcreads


atlas run assembly -w $WD $@

atlas run genomes -w $WD $@
