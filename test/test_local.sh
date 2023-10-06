#! /bin/bash

set -euo pipefail

test_script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
debug_dir="$test_script_dir/../../Debug_atlas"

mkdir -p $debug_dir
cd $debug_dir


reads_dir="test_reads"

snakemake_args=" --quiet rules $@ " 

# if test_reads doean't exist download it

if [ ! -d "$reads_dir" ]; then
    echo "Downloading test reads"
    wget https://zenodo.org/record/3992790/files/test_reads.tar.gz
    tar -xzf test_reads.tar.gz 
    rm test_reads.tar.gz
fi

WD='wd'


rm -f $WD/samples.tsv

#
atlas init $reads_dir -w $WD 

atlas run None screen -w $WD qc  $snakemake_args

echo "\n\nFinished screen\n\n"

atlas run -w $WD qc  $snakemake_args

echo "\n\nFinished qc\n\n"


atlas run assembly -w $WD $snakemake_args

echo "\n\nFinished assembly\n\n"

atlas run binning -w $WD $snakemake_args

echo "\n\nFinished binning\n\n"

atlas run genecatalog --omit-from combine_egg_nogg_annotations combine_dram_genecatalog_annotations -w $WD $snakemake_args

echo "\n\nFinished genecatalog\n\n"

# atlas run genomes -w $WD $@

# echo "\n\nFinished genomes\n\n"

# atlas run strains -w $WD $@

# echo "\n\nFinished strains\n\n"



