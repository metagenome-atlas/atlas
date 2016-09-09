# ATLAS

# Install

```
bedtools
pysam
flash
trimmomatic
megahit
bowtie2
fastqc
maxbin
	hmmer3
	idba-ud
	fraggenescan (1.30 - which actually doesn't exist)
```

# Current Protocol

1. Build contaminant DBs
2. Join R1 and R2
3. Decontaminate across DBs
4. Trim the merged reads
5. FastQC on passing trimmed, R1 failures, and R2 failures
6. Interleave R1 and R2 failures
7. Concatenate #4 output with #6 output
8. Assemble
9. Length filter assembled contigs
10. Statistics (# contigs, total length, GC, etc.)
11. Binning with contigs, interleaved, and non-joined reads

# TODO

- [ ] Gene calling across final contigs
- [ ] Full report on assembly, bins, and genes
- [ ] Incorporate Prokka-like functionality to perform above step and annotate


# Proposed Usage

```
snakemake --jobs 24 \
	--configfile config/atlas.yaml \
	--config eid=test-data
```

Preparing to run, place FASTQ files into `results/<eid>/demultiplexed`. At this point, lets say they're all not compressed and are formatted:

```
<sample-name-1>_R1.fastq
<sample-name-1>_R2.fastq
<sample-name-2>_R1.fastq
<sample-name-2>_R2.fastq
```
