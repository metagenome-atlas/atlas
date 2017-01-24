.. _example-configuration:

Example Configuration
=====================

Lines starting with '#' are treated as comments::

    samples:
        Sample-1:
            path:
                - /data/sample-1_R1.fastq.gz
                - /data/sample-1_R2.fastq.gz
        Sample-2:
            path:
                - /data/sample-2_interleaved.fastq.gz
            type: metatranscriptome
        Sample-3:
            path:
                - /data/sample-3_se.fastq.gz
            paired: false

    threads: 23

    preprocessing:
        # fastas to parse into kmers which are then trimmed from the sequences
        adapters: /pic/projects/mint/atlas_databases/adapters.fa
        # Look for shorter kmers at read tips down to this length; 0 disables
        mink: 8
        minimum_base_quality: 10
        # kmer mismatches allowed during adapter trim process
        allowable_kmer_mismatches: 1
        # length of kmer to search against sequences
        reference_kmer_match_length: 31
        # passing single-end read length
        minimum_passing_read_length: 51
        # Discard reads with a minimum base frequency below this
        min_base_frequency: 0.05

        contamination:
            references:
                # rRNA is a required key here
                rRNA: /pic/projects/mint/atlas_databases/silva_rfam_all_rRNAs.fa
                phiX: /pic/projects/mint/atlas_databases/phiX174_virus.fa
            # Don't look for indels longer than this
            maxindel: 20
            # Fraction of max alignment score required to keep a site
            minratio: 0.65
            # mapping kmer length; range 8-15; longer is faster but uses more memory; shorter is more sensitive
            k: 12
            # Minimum number of seed hits required for candidate sites
            minhits: 1
            # Set behavior on ambiguously-mapped reads (with multiple top-scoring mapping locations):
            #   best    (use the first best site)
            #   toss    (consider unmapped)
            #   random  (select one top-scoring site randomly)
            #   all     (retain all top-scoring sites)
            ambiguous: best

        normalization:
            # kmer length
            k: 21
            # target normalization depth
            t: 100
            # reads must have at least this many kmers over min depth to be retained
            minkmers: 8

    assembly:
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
        spades_k: auto

        # Discard contigs with lower average coverage.
        minc: 5
        # Discard contigs with a lower percent covered bases.
        minp: 40
        # Discard contigs with fewer mapped reads.
        minr: 0
        # Discard contigs shorter than this (after trimming).
        minl: 200
        # Trim the first and last X bases of each sequence.
        trim: 0

    annotation:
        # when counting reads aligning to ORFs, require at least this many bp
        # overlapping the ORF
        minimum_overlap: 20

        references:
            eggnog:
                namemap: /pic/projects/mint/atlas_databases/eggnog.db
                dmnd: /pic/projects/mint/atlas_databases/eggnog.dmnd
                chunk_size: 250000
                # setting top_seqs to 5 will report all alignments whose score is
                # at most 5% lower than the top alignment score for a query
                top_seqs: 5
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
                # filters ORF BLAST hits by only keep hits within this fraction of
                # the highest bitscore; this is recommended over max_hits
                top_fraction: 0.50
            refseq:
                namemap: /pic/projects/mint/atlas_databases/refseq.db
                tree: /pic/projects/mint/atlas_databases/refseq.tree
                dmnd: /pic/projects/mint/atlas_databases/refseq.dmnd
                summary_method: best
                aggregation_method: lca-majority
                majority_threshold: 0.51
                min_length: 60
                max_hits: 10
                top_fraction: 0.50
            expazy:
                namemap: /pic/projects/mint/atlas_databases/expazy.db
                dmnd: /pic/projects/mint/atlas_databases/expazy.dmnd
            cazy:
                namemap: /pic/projects/mint/atlas_databases/cazy.db
                dmnd: /pic/projects/mint/atlas_databases/cazy.dmnd
            cog:
                namemap: /pic/projects/mint/atlas_databases/cog.db
                dmnd: /pic/projects/mint/atlas_databases/cog.dmnd

    summary_counts:
        taxonomy:
            levels:
                - phylum
                - class
                - order
                - species
            COG:
                - cog_id
                - cog_functional_class
                - cog_annotation
            CAZy_EC:
                - cazy_ec
            CAZy_family:
                - cazy_family
            ExPAZy:
                - expazy_name
                - expazy_ec
            RefSeq:
                - refseq_product
        CAZY_EC:
            - cazy_ec
        COG:
            - cog_id
            - cog_functional_class
            - cog_annotation
        ExPAZy:
            - expazy_name
            - expazy_ec
        RefSeq:
            - refseq_product
