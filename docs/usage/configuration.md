\_configuration:

# Configure Atlas

\_contaminants:

## Remove reads from Host

One of the most important steps in the Quality control is to remove
reads from the host\'s genome. You can add any number of genomes to be
removed.

We recommend using genomes where repetitive sequences are masked. See
here for more details [human
genome](http://seqanswers.com/forums/archive/index.php/t-42552.html).

## Co-abundance Binning

::: {#cobinning}
While binning each sample individually is faster, using co-abundance for
binning is recommended. Quantifying the coverage of contigs across
multiple samples provides valuable insights about contig co-variation.
:::

There are two primary strategies for co-abundance binning:

1.  **Cross mapping:** Map the reads from multiple samples to each
    sample\'s contigs.
2.  **Co-binning:** Concatenate contigs from multiple samples and map
    all the reads to these combined contigs.

`final_binner: metabat2` is used for cross-mapping, while
`vamb` or `SemiBin` is used for co-binning.

The samples to be binned together are specified using the
`BinGroup` in the `sample.tsv` file. The size of
the BinGroup should be selected based on the binner and the co-binning
strategy in use.

Cross-mapping complexity scales quadratically with the size of the
BinGroup since each sample's reads are mapped to each other. This might
yield better results for complex metagenomes, although no definitive
benchmark is known. On the other hand, co-binning is more efficient, as
it maps a sample's reads only once to a potentially large assembly.

### Default Behavior

Starting with version 2.18, Atlas places every sample in a single
BinGroup and defaults to `vamb` as the binner unless there
are very few samples. For fewer than 8 samples, `metabat` is
the default binner.

::: note
::: title
Note
:::

This represents a departure from previous versions, where each sample
had its own BinGroup. Running `vamb` in those versions would
consider all samples, regardless of their BinGroup. This change might
cause errors if using a `sample.tsv` file from an older
Atlas version. Typically, you can resolve this by assigning a unique
BinGroup to each sample.
:::

The mapping threshold has been adjusted to 95% identity (single sample
binning is 97%) to allow reads from different strains --- but not other
species --- to map to contigs from a different sample.

If you're co-binning more than 150-200 samples or cross-mapping more
than 50 samples, Atlas will issue a warning regarding excessive samples
in a BinGroup. Although VAMB's official publication suggests it can
handle up to 1000 samples, this demands substantial resources.

Therefore, splitting your samples into multiple BinGroups is
recommended. Ideally, related samples, or those where the same species
are anticipated, should belong to the same BinGroup.

### Single-sample Binning

To employ single-sample binning, simply assign each sample to its own
BinGroup and select `metabat` or `DASTool` as
the `final_binner`.

Although it's not recommended, it's feasible to use
`DASTool` and feed it inputs from `metabat` and
other co-abundance-based binners.

Add the following lines to your \`config.yaml\`:

``` yaml
final_binner: DASTool

binner:
  - metabat
  - maxbin
  - vamb
```

## Long reads {#longreads}

Limitation: Hybrid assembly of long and short reads is supported with
spades and metaSpades. However, metaSpades needs a paired-end short-read
library.

The path of the (preprocessed) long reads should be added manually to
the sample table under a new column heading \'longreads\'.

In addition, the type of the long reads should be defined in the config
file: `longread_type` one of \[\"pacbio\", \"nanopore\", \"sanger\",
\"trusted-contigs\", \"untrusted-contigs\"\]

## Example config file

```yaml
###################################################################
####                 _______   _                    _____      ####
####         /\     |__   __| | |          /\      / ____|     ####
####        /  \       | |    | |         /  \    | (___       ####
####       / /\ \      | |    | |        / /\ \    \___ \      ####
####      / ____ \     | |    | |____   / ____ \   ____) |     ####
####     /_/    \_\    |_|    |______| /_/    \_\ |_____/      ####
####                                                           ####
###################################################################

#  For more details about the config values see:
#  https://metagenome-atlas.rtfd.io

########################
# Execution parameters
########################
# threads and memory (GB) for most jobs especially from BBtools, which are memory demanding
threads: 8
mem: 60

# threads and memory for jobs needing high amount of memory. e.g GTDB-tk,checkm or assembly
large_mem: 250
large_threads: 16
assembly_threads: 8
assembly_memory: 250
simplejob_mem: 10
simplejob_threads: 4

#Runtime only for cluster execution
runtime: #in h
  default: 5
  assembly: 48
  long: 24
  simplejob: 1

# directory where databases are downloaded with 'atlas download'
database_dir: databases

########################
# Quality control
########################
data_type: metagenome # metagenome or metatranscriptome
interleaved_fastqs: false

# remove (PCR)-duplicated reads using clumpify
deduplicate: true
duplicates_only_optical: false
duplicates_allow_substitutions: 2

# used to trim adapters from reads and read ends
preprocess_adapters: /path/to/databases/adapters.fa
preprocess_minimum_base_quality: 10
preprocess_minimum_passing_read_length: 51
# 0.05 requires at least 5 percent of each nucleotide per sequence
preprocess_minimum_base_frequency: 0.05
preprocess_adapter_min_k: 8
preprocess_allowable_kmer_mismatches: 1
preprocess_reference_kmer_match_length: 27
# error correction where PE reads overlap
error_correction_overlapping_pairs: true
#contamination references can be added such that -- key: /path/to/fasta
contaminant_references:
  PhiX: /path/to/databases/phiX174_virus.fa
#  host:/path/to/host_genome.fasta

# We won't allow large indels
contaminant_max_indel: 20
contaminant_min_ratio: 0.65
contaminant_kmer_length: 13
contaminant_minimum_hits: 1
contaminant_ambiguous: best

########################
# Pre-assembly-processing
########################

# Advanced Error correction
error_correction_before_assembly: true
spades_skip_BayesHammer: true # Skip error correction in spades assembler
error_correction_kmer: 31 # can be longer e.g. 62 but takes more memory

# remove reads with k-mers that cannot be used for assembly.
# Filter reads that have a 10% of k-mers below a minimum depth.
error_correction_remove_lowdepth: false
error_correction_minimum_kmer_depth: 1 #
error_correction_aggressive: false

# Merging of pairs
# join R1 and R2 at overlap; unjoined reads are still utilized
merge_pairs_before_assembly: true
merging_k: 62

########################
# Assembly
########################
# megahit OR spades
assembler: spades

minimum_contig_length: 1000
# Megahit
#-----------
# 2 is for metagenomes, 3 for genomes with 30x coverage
megahit_min_count: 2
megahit_k_min: 21
megahit_k_max: 121
megahit_k_step: 20
megahit_merge_level: 20,0.98
megahit_prune_level: 2
megahit_low_local_ratio: 0.2
# ['default','meta-large','meta-sensitive']
megahit_preset: default

# Spades
#------------
spades_use_scaffolds: true # if false use contigs
#Comma-separated list of k-mer sizes to be used (all values must be odd, less than 128 and listed in ascending order).
spades_k: auto
spades_preset: meta # meta, ,normal, rna  single end libraries doesn't work for metaspades
spades_extra: ""
longread_type: none # [none,"pacbio", "nanopore", "sanger", "trusted-contigs", "untrusted-contigs"]
# Preprocessed long reads can be defined in the sample table with 'longreads' , for more info see the spades manual

# Filtering
#------------
# filter out assembled noise
# this is more important for assembly from megahit
filter_contigs: false
# trim contig tips
contig_trim_bp: 0
# require contigs to have read support
minimum_average_coverage: 1
minimum_percent_covered_bases: 20
minimum_mapped_reads: 0

########################
# Quantification
########################

# Mapping reads to contigs
#--------------------------
contig_min_id: 0.9
contig_map_paired_only: true
contig_max_distance_between_pairs: 1000
maximum_counted_map_sites: 10
minimum_map_quality: 0

########################
# Binning
########################

final_binner: vamb # [SemiBin, vamb, metabat, DASTool]

semibin_options: ""

metabat:
  sensitivity: sensitive
  min_contig_length: 1500 # metabat needs >1500

maxbin:
  max_iteration: 50
  prob_threshold: 0.9
  min_contig_length: 1000

DASTool:
  search_engine: "diamond"
  score_threshold: 0.5 # Score threshold until selection algorithm will keep selecting bins [0..1].

genome_filter_criteria: "(Completeness-5*Contamination >50 ) & (Length_scaffolds >=50000) & (Ambigious_bases <1e6) & (N50 > 5*1e3) & (N_scaffolds < 1e3)"

filter_chimieric_bins: true # filter chimeric bins using GUNC
gunc_database: "progenomes" # 'progenomes' or 'gtdb'

genome_dereplication:
  ANI: 0.95 ## Genome dreplication threshold 0.95 is more or less species
  overlap: 0.2

rename_mags_contigs: true #Rename contigs of representative MAGs

########################
# Annotations
#######################

annotations:
  - gtdb_tree
  - gtdb_taxonomy
  - genes
  - kegg_modules
  - dram

########################
# Gene catalog
#######################
genecatalog:
  source: contigs # [contigs, genomes] Predict genes from all contigs or only from the representative genomes
  clustermethod: linclust # [mmseqs or linclust] see mmseqs for more details
  minlength_nt: 270 # min length
  minid: 0.90 # min id for gene clustering for the main gene catalog used for annotation
  coverage: 0.9
  extra: " "
  SubsetSize: 500000

gene_annotations:
  - eggNOG
  # - dram

eggNOG_use_virtual_disk: false # coping the eggNOG DB to a virtual disk can speed up the annotation
virtual_disk: "/dev/shm" # But you need 37G extra ram
```

## Detailed configuration

toctree::

:

    maxdepth

    :   1

    ../advanced/qc ../advanced/assembly
