


# Change log


## 2.18.1

Fix error with downloading DRAM. Update to DRAM v1.5

## 2.18

- Qc reads, assembly are now written in the sample.tsv from the start. This should fix errors of partial writing to the sample.tsv https://github.com/metagenome-atlas/atlas/issues/695
- It also allows you to add external assemblies.
- singletons reads are no longer used trough the pipeline. 
- This changes the default paths for raw reads and assemblies. 
assembly are now in `Assembly/fasta/{sample}.fasta`
reads: `QC/reads/{sample}_{fraction}.fastq.gz`

**Seemless update**: If you update atlas and continue on an old project. Your old files will be copied.
Or the path defined in the sample.tsv will be used. 



## 2.17

### Skani
The tool Skani claims to be better and faster than the combination of mash + FastANI as used by dRep
I implemented the skin for species clustering. 
We now do the species clustering in the `atlas run binning` step. 
So you get information about the number of dereplicated species in the binning report. This allows you to run different binners before choosing the one to use for the genome annotation. 
Also, the file storage was improved all important files are in `Binning/{binner}/`



My custom species clustering does the following steps:

1. Pre-cluster genomes with *single-linkage* at 92.5 ANI. 
2. **Re-calibrate checkm2 results.**
 - If a minority of genomes from a pre-cluster use a different translation table they are removed
 - If some genomes of a pre-cluster don't use the specialed completeness model we re-calibrate completeness to the minimum value.
 This ensures that not a bad genome evaluated on the general model is preferred over a better genome evaluated on the specific model.
See also https://silask.github.io/post/better_genomes/ Section 2.
- Drop genomes that don't correspond to the filter criteria after re-calibration
3. Cluster genomes with ANI threshold default 95%
4. Select the best genome as representative based on the Quality score Completeness - 5x Contamination





### New Contributors
* @jotech made their first contribution in https://github.com/metagenome-atlas/atlas/pull/667

## 2.16

 * gtdb08

## 2.15

* Use Gunc
* New Folder organisation: Main output files for Binning are in the new folder `Binning`
* Use hdf-format for gene catalogs. Allow efficient storage and selective access to large count and coverage matrices from the genecatalog. (See docs for how to load them) https://github.com/metagenome-atlas/atlas/pull/621
* Semibin v. 1.5 by @SilasK in https://github.com/metagenome-atlas/atlas/pull/622


## 2.14

* Support for checkm2  by @SilasK in https://github.com/metagenome-atlas/atlas/pull/607

Thank you @trickovicmatija for your help.

**Full Changelog**: https://github.com/metagenome-atlas/atlas/compare/v2.13.1...v2.14.0
## 2.13

* use minimap for contigs, genecatalog and genomes in https://github.com/metagenome-atlas/atlas/pull/569  https://github.com/metagenome-atlas/atlas/pull/577
* filter genomes my self  in https://github.com/metagenome-atlas/atlas/pull/568
The filter function is defined in the config file:
```
genome_filter_criteria: "(Completeness-5*Contamination >50 ) & (Length_scaffolds >=50000) & (Ambigious_bases <1e6) & (N50 > 5*1e3) & (N_scaffolds < 1e3)"
```
The genome filtering is similar as other publications in the field, e.g. GTDB. What is maybe a bit different is that genomes with completeness around 50% **and** contamination around 10% are excluded where as using the default parameters dRep would include those. 

* use Drep again  in https://github.com/metagenome-atlas/atlas/pull/579
We saw better performances using drep. This scales also now to ~1K samples
* Use new Dram version 1.4 by in https://github.com/metagenome-atlas/atlas/pull/564


**Full Changelog**: https://github.com/metagenome-atlas/atlas/compare/v2.12.0...v2.13.0

## 2.12

* GTDB-tk requires rule `extract_gtdb` to run first by @Waschina in https://github.com/metagenome-atlas/atlas/pull/551
* use Galah instead of Drep 
* use bbsplit for mapping to genomes (maybe move to minimap in future)
* faster gene catalogs quantification using minimap. 
* Compatible with snakemake v7.15
### New Contributors
* @Waschina made their first contribution in https://github.com/metagenome-atlas/atlas/pull/551

**Full Changelog**: https://github.com/metagenome-atlas/atlas/compare/v2.11.1...v2.12.0

## 2.11
* Make atlas handle large gene catalogs using parquet and pyfastx (Fix #515)

parquet files can be opened in python with 
```
import pandas as pd
coverage = pd.read_parquet("working_dir/Genecatalog/counts/median_coverage.parquet")
coverage.set_index("GeneNr", inplace=True)

```

and in R it should be something like:

```
arrow::read_parquet("working_dir/Genecatalog/counts/median_coverage.parquet")

```


**Full Changelog**: https://github.com/metagenome-atlas/atlas/compare/v2.10.0...v2.11.0

## [2.10](https://github.com/metagenome-atlas/atlas/compare/v2.9.1...v2.10.0) 

### Features
* GTDB version 207
* Low memory taxonomic annotation


## [2.9](https://github.com/metagenome-atlas/atlas/compare/v2.8.2...v2.9.0) 

### Features
* ✨ Start an atlas project from public data in SRA [Docs](https://metagenome-atlas.readthedocs.io/en/latest/usage/getting_started.html#start-a-new-project-with-public-data)
* Make atlas ready for python 3.10  https://github.com/metagenome-atlas/atlas/pull/498
* Add strain profiling using inStrain You can run `atlas run genomes strains`

### New Contributors
* @alienzj made their first contribution to fix config when run DRAM annotate in https://github.com/metagenome-atlas/atlas/pull/495


## 2.8
This is a major update of metagenome-atlas. It was developed for the [3-day course in Finnland](https://silask.github.io/talk/3-day-course-on-metagenome-atlas/), that's also why it has a finish release name. 


### New binners
It integrates bleeding-edge binners `Vamb` and `SemiBin` that use Co-binning based on co-abundance. Thank you @yanhui09 and @psj1997 for helping with this. The first results show better results using these binners over the default. 

[See more](https://metagenome-atlas.readthedocs.io/en/v2.8.0/usage/output.html#binning)

### Pathway annotations
The command `atlas run genomes` produces genome-level functional annotation and Kegg pathways respective modules. It uses DRAM  from @shafferm with a hack to produce all available Kegg modules. 

[See more](https://metagenome-atlas.readthedocs.io/en/v2.8.0/usage/output.html#annotations)

### Genecatalog
The command `atlas run genecatalog` now produces directly the abundance of the different genes. See more in #276 

> In future this part of the pipeline will include protein assembly to better tackle complicated metagenomes. 

### Minor updates

#### Reports are back
See for example the [QC report](https://metagenome-atlas.readthedocs.io/en/v2.8.0/_static/QC_report.html)

#### Update of all underlying tools
All tools use in atlas are now up to date.  From assebler to GTDB.
The one exception is, BBmap which contains a [bug](https://sourceforge.net/p/bbmap/tickets/48/) and ignores the minidenty parameter.

#### Atlas init 
Atlas init correctly parses fastq files even if they are in subfolders and if paired-ends are named simply  Sample_1/Sample_2. @Sofie8 will be happy about this.
Atlas log uses nice colors.
 
#### Default clustering of Subspecies

The default ANI threshold for genome-dereplication was set to 97.5% to include more sub-species diversity. 

[See more](https://metagenome-atlas.readthedocs.io/en/v2.8.0/usage/output.html#genomes) 








