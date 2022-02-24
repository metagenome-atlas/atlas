#! /bin/bash
set -euo pipefail

atlas --version

Test_dir="test/Test_sra_init"

rm -rf $Test_dir
mkdir -p $Test_dir

echo "Download reads from our library"

WD=$Test_dir/"Mouse"

atlas init-public PRJEB20796 -w $WD --overwrite

echo "Run Atlas"

atlas run qc -w $WD --dry-run $@


echo "Download reads from HMP"
WD=$Test_dir/"HMP"

# this fails as HMP have samples sequenced with different platforms
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

 


echo "Run Atlas"
atlas init-public SRP002423 -w $WD


atlas run qc -w $WD --dry-run $@



## smal data


echo "Download reads from small dataset for real test"

WD=$Test_dir/"Small"

echo "gives error as library is selected with PCR"
set +e
atlas init-public SAMEA9831203 SAMEA9831204 -w $WD
set -e

echo "use this data anyway"
cp $WD/SRA/RunInfo_original.tsv $WD/RunInfo.tsv

atlas init-public Anything -w $WD 

echo "Run Atlas"

atlas run None download_sra -w $WD $@