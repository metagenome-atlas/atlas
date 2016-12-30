<p align="center">
    <img src=resources/images/logo.jpg width=300 />
</p>

# ATLAS

# Install

All dependencies are installed via [`conda`](https://www.continuum.io/downloads) using the [bioconda](https://github.com/bioconda/bioconda-recipes) channel. The workflow and some dependencies require Python 3.5.

To install:

```
conda install -c bioconda \
    bbmap diamond fastqc megahit prodigal samtools snakemake spades verse
```

Then, as `atlas` is still in active development:

```
git clone git@github.com:????????/atlas.git
```

## Databases

To download the databases and their respective metadata databases...

# Usage

```
snakemake --configfile config.yaml
```

# Configuration

## Samples

Samples are defined with a name, file path(s) and the type of data. A single file path is interpreted as interleaved paired-end reads, while two paths must include full paths to R1 and R2.

For `type`, the value can be either `metagenome` or `metatranscriptome`. The default if not specified is `metagenome`.

Interleaved PE input:

```
samples:
    sample-1:
        path:
            - /data/sample-1_pe.fastq.gz
        type: metagenome
```

Paired-end as separate files:

```
samples:
    sample-1:
        path:
            - /data/sample-1_R1.fastq.gz
            - /data/sample-1_R2.fastq.gz
        type: metagenome
```

## Temporary Directory

Some steps, like assembly, may have an optional temporary directory. If specified in the configuration or as an environmental variable [TMPDIR|TEMP|TMP], the step will use this directory rather than the current working directory.

**Default: the working directory**

```
temporary_directory: /scratch
```

## Threads

Most steps of the workflow are utilizing applications that can thread or otherwise use multiple cores. Leaving this one below the max, in cases where many samples are being analyzed, may be optimal as single-threaded jobs will be processed more efficiently.

Threads is used per step and will most likely be a subset of `--jobs`. `--jobs` represents the total available cores for all simultaneous steps.

**Default: 1**

```
threads: 23
```

## Shell Command Prefix

This allows the user to prepend something to each command. In our case, we reserve a `slurm` allocation, then submit jobs across the reservation with:

```
prefix: srun --exclusive -N1 -n1 -c
```

And `threads` is appended to the prefix when the workflow is executed.

**Default: None**

## Preprocessing of Reads

This starts the preprocessing subsection of the configuration file. All settings will be indented from "preprocessing" and indentation levels matter.

```
preprocessing:
    adapters: /databases/adapters.fa
    minimum_base_quality: 10
    min_base_frequency: 0.05
    contamination:
        references:
            rRNA: /refs/rrna.fasta
        k: 12
        ambiguous: best
    normalization:
        k: 21
        t: 100
```

### Adapters

FASTA file paths for adapter sequences to be trimmed from the sequence ends. It is best practice to use full file paths to reference files.

We provide the adapter reference FASTA included in `bbmap`.

```
preprocessing:
    adapters: /databases/adapters.fa
```

### Quality Trimming

Trim regions with an average quality below this threshold. Higher is more stringent.

**Default: 10**

```
preprocessing:
    minimum_base_quality: 10
```

### Adapter Trimming at Read Tips

Allow shorter kmer matches down to `mink` at the read ends. 0 disables.

**Default: 8**

```
preprocessing:
    mink: 8
```

### Allowable Mismatches in Adapter Hits

Maximum number of substitutions between the target adapter kmer and the query sequence kmer. Lower is more stringent.

**Default: 1**

```
preprocessing:
    allowable_kmer_mismatches: 1
```

### Kmer Length

Kmer length used for finding contaminants. Contaminant matches shorter than this length will not be found.

**Default: 27**

```
preprocessing:
    reference_kmer_match_length: 27
```

### Read Length Threshold

This is applied after quality and adapter trimming have been applied to the sequence.

**Default: 51**

```
preprocessing:
    minimum_passing_read_length: 51
```

### Complexity Filter

Require this fraction of each nucleotide per sequence to eliminate low complexity reads.

**Default: 0.05**

```
preprocessing:
    min_base_frequency: 0.05
```

### Contamination Parameters

Contamination reference sequences in the form of nucleotide FASTA files can be provided and filtered from the reads using the following parameters. These still fall within the 'preprocessing' section of the configuration.

#### Maximum Insertion/Deletion

Have `bbsplit.sh` stop searching for possible mappings with indels longer than this. Lower is faster.

**Default: 20**

```
preprocessing:
    contamination:
        maxindel: 20
```

#### Required Mapped Read Fraction

Of the possible maximum alignment score, force at least this fraction per mapping.

**Default: 0.65**

```
preprocessing:
    contamination:
        minratio: 0.65
```

#### Minimum Seed Hits

Minimum number of seed hits required for candidate sites.

**Default: 1**

```
preprocessing:
    contamination:
        minhits: 1
```

#### Ambiguous Mappings

The method for which we will deal with reads that map to multiple contamination reference sequences. Possible values include:

| Value  | Definition                            |
|--------|---------------------------------------|
| best   | Use the first best site.              |
| toss   | Consider the read unmapped.           |
| random | Select one top-scoring site randomly. |
| all    | Retain all top-scoring sites.         |

**Default: best**

```
preprocessing:
    contamination:
        ambiguous: best
```

#### Mapping Kmer Length

Mapping kmer length in the range of 8 to 15. Shorter will be more sensitive and slower.

**Default: 13**

```
preprocessing:
    contamination:
        k: 13
```

#### Reference Sequences

##### Ribosomal RNA (`rRNA`)

This reference FASTA is required though you can provide an alternate to the provided rRNA reference.

##### Additional References

Any number of additional contamination reference sequences can be used. The key is the name that will be integrated into the file name and provide the path to the file such that:

```
preprocessing:
    contamination:
        references:
            rRNA: /refs/rrna.fasta
            human: /refs/human.fasta
            cat: /refs/cat.fasta
```

### Normalization Parameters

#### Kmer Length (`k`)

#### Target Coverage (`t`)

#### Minimum Passing Kmers (`minkmers`)

Reads must have at least this many kmers over the minimum depth to be retained.

**Default: 8**

assembly:
    # 'spades' or 'megahit'
    assembler: megahit
    # fraction of the machine's total memory or bytes
    memory: 0.99
    # minimum multiplicity for filtering (k_min+1)-mers
    minimum_count: 2
    # minimum kmer size (<= 255), must be odd number
    kmer_min: 21
    # maximum kmer size (<= 255), must be odd number
    kmer_max: 121
    # increment of kmer size of each iteration (<= 28), must be even number
    kmer_step: 20
    # merge complex bubbles of length <= l*kmer_size and similarity >= s
    merge_level: 20,0.98
    # strength of low depth pruning (0-3)
    prune_level: 2
    # ratio threshold to define low local coverage contigs
    low_local_ratio: 0.2
    # minimum length of contigs to output from the assembler; can be filtered
    # downstream using minl
    minimum_contig_length: 200
    # comma-separated list of k-mer sizes (must be odd and less than 128)
    spades_k: auto
    # Discard contigs with lower average coverage.
    minc: 5
    # Discard contigs with a lower percent covered bases.
    minp: 40
    # Discard contigs with fewer mapped reads.
    minr: 0
    # Discard contigs shorter than this (after trimming).
    minl: 250
    # Trim the first and last X bases of each sequence.
    trim: 0

annotation:
    ## ORFs
    # https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    translation_table: 11
    # when counting reads aligning to ORFs, require at least this many bp
    # overlapping the ORF
    minimum_overlap: 20

    references:
        eggnog:
            # non-tree based reference requires namemap database and fasta
            namemap: /pic/projects/mint/atlas_databases/functional/eggnog/eggnog4_nonredundant.db
            fasta: /pic/projects/mint/atlas_databases/functional/eggnog/eggnog4_nonredundant.fasta
            # number of entries per FASTA to be aligned with DIAMOND
            chunk_size: 250000
            # 'fast' or 'sensitive'
            run_mode: fast
            # setting top_seqs to 5 will report all alignments whose score is
            # at most 5% lower than the top alignment score for a query
            top_seqs: 5
            # maximum e-value to report alignments
            e_value: "0.000001"
            # minimum identity % to report an alignment
            min_identity: 50
            # minimum query cover % to report an alignment
            query_coverage: 60
            # gap open penalty
            gap_open: 11
            # gap extension penalty
            gap_extend: 1
            # Block size in billions of sequence letters to be processed at a time.
            # This is the main parameter for controlling DIAMOND's memory usage.
            # Bigger numbers will increase the use of memory and temporary disk space,
            # but also improve performance. The program can be expected to roughly use
            # six times this number of memory (in GB).
            block_size: 4
            # The number of chunks for processing the seed index (default=4). This
            # option can be additionally used to tune the performance. It is
            # recommended to set this to 1 on a high memory server, which will
            # increase performance and memory usage, but not the usage of temporary
            # disk space.
            index_chunks: 4
            # 'majority' or 'best'; summary method for annotating ORFs
            summary_method: best
            # minimum allowable BLAST alignment length
            min_length: 60
            # maximum allowable e-value of BLAST hit when parsing DIAMOND hits
            max_evalue: 0.000001
            # maximum number of BLAST hits to consider when summarizing ORFs
            max_hits: 10
            # filters ORF BLAST hits by only keep hits within this fraction of
            # the highest bitscore; this is recommended over max_hits
            top_fraction: 0.50
            # minimum allowable BLAST alignment bitscore; 0 effectively disables
            min_bitscore: 0
        refseq:
            # tree based reference requires namemap database, tree, and fasta
            namemap: /pic/projects/mint/atlas_databases/functional/refseq/refseq78.complete.nonredundant_protein.faa.db
            tree: /pic/projects/mint/atlas_databases/functional/refseq/refseq78.complete.nonredundant_protein.faa.tree
            fasta: /pic/projects/mint/atlas_databases/functional/refseq/refseq78.complete.nonredundant_protein.faa
            # number of entries per FASTA to be aligned with DIAMOND
            chunk_size: 250000
            run_mode: fast
            top_seqs: 5
            e_value: "0.000001"
            min_identity: 50
            query_coverage: 60
            gap_open: 11
            gap_extend: 1
            block_size: 6
            index_chunks: 1
            # 'lca', 'majority', or 'best'; summary method for annotating ORFs; when
            # using LCA, it's recommended that one limits the number of hits using a
            # low top_fraction
            summary_method: best
            # 'lca', 'lca-majority', or 'majority'; summary method for aggregating ORF
            # taxonomic assignments to contig level assignment; 'lca' will result in
            # most stringent, least specific assignments
            aggregation_method: lca-majority
            # constitutes a majority fraction at tree node for 'lca-majority' ORF
            # aggregation method
            majority_threshold: 0.51
            # minimum allowable BLAST alignment length
            min_length: 60
            # maximum allowable e-value of BLAST hit
            max_evalue: 0.000001
            # maximum number of BLAST hits to consider when summarizing ORFs; can
            # drastically alter ORF LCA assignments if too high without further limits
            max_hits: 10
            top_fraction: 0.50
        expazy:
            namemap: /pic/projects/mint/atlas_databases/functional/expazy/expazy.db
            fasta: /pic/projects/mint/atlas_databases/functional/expazy/expazy.fasta
            chunk_size: 500000
            # 'fast' or 'sensitive'
            run_mode: fast
            top_seqs: 2
            index_chunks: 1
            summary_method: majority
        cazy:
            namemap: /pic/projects/mint/atlas_databases/functional/dbcan/dbcan.db
            fasta: /pic/projects/mint/atlas_databases/functional/dbcan/dbcan.fasta
            chunk_size: 500000
            # 'fast' or 'sensitive'
            run_mode: fast
            top_seqs: 2
            index_chunks: 1
            summary_method: majority

summary_counts:
    # Possible columns table column values upon which to aggregate:
        # contig, orf

        # from refseq:
        # taxonomy, orf_taxonomy, refseq_product

        # from eggnog:
        # uniprot_ac, eggnog_ssid_b, eggnog_species_id, uniprot_id, cog_func_id, cog_id,
        # cog_product, cog_level1_code, cog_level1_name, cog_level2_name,
        # ko_id, ko_level1_name, ko_level2_name, ko_level3_id,
        # ko_level3_name, ko_gene_symbol, ko_product, ko_ec

        # from expazy:
        # expazy_name, expazy_ec

        # from cazy (dbcan):
        # cazy_gene, cazy_family, cazy_class, cazy_ec

    # this is a special case to allow for taxon level specification
    taxonomy:
        # limit taxonomy in classification to the depth specified
        # possible values: kingdom, domain, phylum, class, order, family, genus, species
        # all levels if omitted
        levels:
            - phylum
            - class
            - order
            - species
        # tables to generate at these taxonomic levels
        KO:
            - ko_id
            - ko_ec
        COG:
            - cog_id
        CAZy_EC:
            - cazy_ec
        CAZy_family:
            - cazy_family
        ExPAZy:
            - expazy_name
            - expazy_ec
    KO:
        - ko_id
        - ko_gene_symbol
        - ko_product
        - ko_ec
    KO_lvl1:
        - ko_level1_name
    KO_lvl2:
        - ko_level2_name
    KO_lvl3:
        - ko_level3_name
    CAZY_EC:
        - cazy_ec
    COG:
        - cog_id
        - cog_product
    COG_lvl1:
        - cog_level1_name
        - cog_level2_name
    ExPAZy:
        - expazy_name
        - expazy_ec


# Output
