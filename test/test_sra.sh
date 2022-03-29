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


echo "Download reads from HMP"
WD=$Test_dir/"HMP"
echo "WD="$WD

# this fails as HMP have samples sequenced with different platforms

echo "(expect errors)"
set +e
atlas init-public SRP002423 -w $WD

set -e

echo "drop illumina samples"
sed -i.bak '/ILLUMINA/d' $WD/RunInfo.tsv

# modify assembler as spades cannot handle single end reads

# python << END
# from ruamel.yaml import YAML
# yaml = YAML()
# config_file="$WD/config.yaml"
# config= yaml.load(open(config_file))
# config['assembler'] = 'megahit'
# yaml.dump(config, open(config_file, 'w'))
# END


echo "create sample table"
atlas init-public continue -w $WD

echo "Run Atlas"

atlas run qc -w $WD --dry-run $@

## single end

echo "Now with a single end sample"

WD=$Test_dir/"SingleEnd"
echo "WD="$WD

atlas init-public SAMEA104416160  -w $WD

atlas run None download_sra -w $WD $@

## smal data


echo "Download reads from small dataset for real test"

WD=$Test_dir/"Small"
echo "WD="$WD

echo "gives warning as library is selected with PCR"

atlas init-public SAMEA9831203 SAMEA9831204 -w $WD

echo "Run Atlas"

atlas run None download_sra -w $WD $@ 