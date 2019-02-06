#! /bin/bash
set -euo pipefail





atlas --version




databaseDir="databases"
WD='WD'

rm -fr $WD
#
atlas init --db-dir $databaseDir --threads 3  -w $WD ../example_data


atlas run -w $WD qc $@

atlas run assembly -w $WD $@
