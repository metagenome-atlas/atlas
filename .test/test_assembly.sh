#! /bin/bash
set -euo pipefail





atlas --version



samplenames="Mycoplasma Streptococcus"
databaseDir=".test/databases"
WD='.test/Test_assembly'
reads_dir=".test/reads/stub"

ressource_args=" --config simplejob_mem=4 java_mem=4 assembly_mem=4"


# gen randomreads
#very low number only for assembly
#snakemake -s atlas/rules/testing.smk -d $reads_dir --config reads=500


rm -f $WD/samples.tsv
#
atlas init --db-dir $databaseDir --threads 4  -w $WD $reads_dir


atlas run -w $WD qc $ressource_args $@

atlas run assembly -w $WD $ressource_args $@

atlas run assembly -w $WD $@

atlas run assembly -w $WD $@

echo "copy qc reads and assemble"

WD2='.test/Test_assembly_skipQC'
reads_dir=$WD2/"reads"

rm -f $WD2/samples.tsv
mkdir -p $reads_dir
cp $WD/*/sequence_quality_control/*_QC_R?.fastq.gz $reads_dir

atlas init --db-dir $databaseDir --threads 4 --assembler megahit --skip-qc -w $WD2 $reads_dir

atlas run -w $WD2 assembly $ressource_args $@


echo "start from interleaved QC reads"

WD3='.test/Test_assembly_interleved'
reads_dir=$WD3/"reads"

rm -f $WD3/samples.tsv
mkdir -p $reads_dir

for sample in $samplenames ;
do
reformat.sh in=$WD/$sample/sequence_quality_control/${sample}_QC_R1.fastq.gz \
  in2=$WD/$sample/sequence_quality_control/${sample}_QC_R2.fastq.gz out=$reads_dir/${sample}.fastq.gz overwrite=true
done

atlas init --db-dir $databaseDir --threads 4 --skip-qc -w $WD3 $reads_dir

atlas run -w $WD3 assembly $ressource_args interleaved_fastqs=True $@
