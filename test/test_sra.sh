#! /bin/bash
set -euo pipefail

atlas --version

Test_dir="test/Test_sra_init"

rm -rf $Test_dir
mkdir -p $Test_dir

echo "Download reads from our library"

WD=$Test_dir/"Mouse"
echo "WD="$WD

atlas init-public PRJEB20796 -w $WD

echo "Run Atlas"

atlas run qc -w $WD --dry-run $@




## single end

echo "Now with a single end sample"

WD=$Test_dir/"SingleEnd"
echo "WD="$WD

atlas init-public ERR2213683  -w $WD

atlas run qc -w $WD --dry-run $@





# this fails as HMP have samples sequenced with different platforms
# HMP is also single end

echo "Download reads from HMP"
WD=$Test_dir/"HMP"
echo "WD="$WD


set +e
atlas init-public SRP002423 -w $WD

set -e
echo "(expected errors)"


echo "drop illumina samples"
sed -i.bak '/ILLUMINA/d' $WD/RunInfo.csv



echo "Continue public init"
atlas init-public continue -w $WD

echo "Run Atlas"

atlas run qc -w $WD --dry-run $@


## small data


echo "Download reads from small dataset for real test"

WD=$Test_dir/"Small"
echo "WD="$WD

echo "gives warning as library is selected with PCR"

atlas init-public ERR1739691 ERR1739692 -w $WD

echo "Run Atlas"

atlas run None download_sra -w $WD $@ 


# Test single end

echo "Run single end test"

WD=$Test_dir/"SingleEnd"
echo "WD="$WD

atlas run None download_sra -w $WD $@
