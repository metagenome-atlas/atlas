.. _example-configuration:

Example Configuration
=====================

Lines starting with '#' are treated as comments::

    samples:
        # sample ID
        Sample-1:
            fastq:
                - /data/sample-1_R1.fastq.gz
                - /data/sample-1_R2.fastq.gz
        Sample-2:
            fastq: /data/sample-2_interleaved.fastq.gz
            type: metatranscriptome
        Sample-3:
            fastq: /data/sample-3_se.fastq.gz
            type: metagenome
            paired: false

    tmpdir: /scratch
    threads: 24

    # sets threads for all methods that are not single-threaded
    threads: 24
    java_mem: 16g
    # shell command prefix; useful for unconventional job queue situations
    # variables can be replaced from any top-level key in the configuration, e.g.
    # __threads__ is replaced with threads defined above
    # prefix: "srun --exclusive -N1 -n1 -c__threads__"

    # fastas to parse into kmers which are then trimmed from the sequences
    preprocess_adapters: /pic/projects/mint/atlas_databases/adapters.fa
    # Look for shorter kmers at read tips down to this length; 0 disables
    preprocess_adapter_min_k: 8
    preprocess_minimum_base_quality: 10
    # kmer mismatches allowed during adapter trim process
    preprocess_allowable_kmer_mismatches: 1
    # length of kmer to search against sequences
    preprocess_reference_kmer_match_length: 23
    # passing single-end read length, prior to merging
    preprocess_minimum_passing_read_length: 51
    # Discard reads with a minimum base frequency below this
    preprocess_minimum_base_frequency: 0.05
    perform_error_correction: true

    # CONTAMINATION FILTERS -- if not defined, this step is skipped
    contaminant_references:
        # this key ('rRNA') is required
        rRNA: /pic/projects/mint/atlas_databases/silva_rfam_all_rRNAs.fa
        # any number of arbitrary name:file path references may be added
        phiX: /pic/projects/mint/atlas_databases/phiX174_virus.fa
    # Don't look for indels longer than this
    contaminant_max_indel: 20
    # Fraction of max alignment score required to keep a site
    contaminant_min_ratio: 0.65
    # mapping kmer length; range 8-15; longer is faster but uses more memory; shorter is more sensitive
    contaminant_kmer_length: 12
    # Minimum number of seed hits required for candidate sites
    contaminant_minimum_hits: 1
    # Set behavior on ambiguously-mapped reads (with multiple top-scoring mapping locations):
    #   best    (use the first best site)
    #   toss    (consider unmapped)
    #   random  (select one top-scoring site randomly)
    #   all     (retain all top-scoring sites)
    contaminant_ambiguous: best

    normalization_kmer_length: 21
    normalization_target_depth: 100
    # reads must have at least this many kmers over min depth to be retained
    normalization_minimum_kmers: 8

    # ASSEMBLY
    # 'spades' or 'megahit'
    assembler: megahit
    # fraction of the machine's total memory or bytes
    megahit_memory: 0.99
    # minimum multiplicity for filtering (k_min+1)-mers
    megahit_min_count: 2
    # minimum kmer size (<= 255), must be odd number
    megahit_k_min: 21
    # maximum kmer size (<= 255), must be odd number
    megahit_k_max: 121
    # increment of kmer size of each iteration (<= 28), must be even number
    megahit_k_step: 20
    # merge complex bubbles of length <= l*kmer_size and similarity >= s
    megahit_merge_level: 20,0.98
    # strength of low depth pruning (0-3)
    megahit_prune_level: 2
    # ratio threshold to define low local coverage contigs
    megahit_low_local_ratio: 0.2
    # minimum length of contigs (after contig trimming)
    minimum_contig_length: 200
    # comma-separated list of k-mer sizes (must be odd and less than 128)
    spades_k: auto
    # Discard contigs with lower average coverage.
    minimum_average_coverage: 5
    # Discard contigs with a lower percent covered bases.
    minimum_percent_covered_bases: 40
    # Discard contigs with fewer mapped reads.
    minimum_mapped_reads: 0
    # Trim the first and last X bases of each sequence.
    contig_trim_bp: 0


    ################
    ## ANNOTATION ##
    ################

    # https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    translation_table: 11
    # when counting reads aligning to ORFs, require at least this many bp
    # overlapping the ORF
    minimum_region_overlap: 1
    # count primary read alignments only when getting counts across ORFs
    primary_only: false
    # allow counting of multi-mapped reads when getting counts across ORFs
    count_multi_mapped_reads: true
    maximum_counted_map_sites: 10

    perform_genome_binning: true
    maxbin_max_iteration: 50
    maxbin_min_contig_length: 500
    maxbin_prob_threshold: 0.9

    refseq_namemap: /pic/projects/mint/atlas_databases/refseq.db
    refseq_tree: /pic/projects/mint/atlas_databases/refseq.tree
    diamond_db: /pic/projects/mint/atlas_databases/refseq.dmnd
    # 'fast' or 'sensitive'
    diamond_run_mode: fast
    # setting top_seqs to 5 will report all alignments whose score is
    # at most 5% lower than the top alignment score for a query
    diamond_top_seqs: 2
    # maximum e-value to report alignments
    diamond_e_value: "0.000001"
    # minimum identity % to report an alignment
    diamond_min_identity: 50
    # minimum query cover % to report an alignment
    diamond_query_coverage: 60
    # gap open penalty
    diamond_gap_open: 11
    # gap extension penalty
    diamond_gap_extend: 1
    # Block size in billions of sequence letters to be processed at a time.
    # This is the main parameter for controlling DIAMOND's memory usage.
    # Bigger numbers will increase the use of memory and temporary disk space,
    # but also improve performance. The program can be expected to roughly use
    # six times this number of memory (in GB).
    diamond_block_size: 2
    # The number of chunks for processing the seed index (default=4). This
    # option can be additionally used to tune the performance. It is
    # recommended to set this to 1 on a high memory server, which will
    # increase performance and memory usage, but not the usage of temporary
    # disk space.
    diamond_index_chunks: 1
    # 'lca', 'majority', or 'best'; summary method for annotating ORFs; when
    # using LCA, it's recommended that one limits the number of hits using a
    # low top_fraction
    summary_method: lca
    # 'lca', 'lca-majority', or 'majority'; summary method for aggregating ORF
    # taxonomic assignments to contig level assignment; 'lca' will result in
    # most stringent, least specific assignments
    aggregation_method: lca-majority
    # constitutes a majority fraction at tree node for 'lca-majority' ORF
    # aggregation method
    majority_threshold: 0.51
