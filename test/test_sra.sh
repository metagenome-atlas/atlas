#! /bin/bash
set -euo pipefail

atlas --version

Test_dir="test/Test_sra_init"

rm -rf $Test_dir
mkdir -p $Test_dir

echo "Download reads from our library"

WD=$Test_dir/"Mouse"

atlas init-public PRJEB20796 -w $WD



echo "Download reads from HMP"
WD=$Test_dir/"HMP"

# this fails as HMP have samples sequenced with different platforms
set +e
atlas init-public SRP002423 -w $WD

set -e

echo "drop illumina samples"
sed -i.bak '/ILLUMINA/d' $WD/SRA/Runtable.tsv


