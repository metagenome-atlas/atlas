Annotation
==========

Any column later specified in summary tables requires that its reference be
specified. For instance, specifying taxonomic breakdowns across EC will require
RefSeq for the taxonomic assignments and ENZYME for the EC data.

Minimal working configuration:

featureCounts Mode:

-z 0   featureCounts(default): Quantify by overlapping and voting.
           If the read pair overlaps multiple genes, it will assign the read pair
           to the gene that is overlapped by both reads.

           Please refer to the SUBREAD users guide for more information.
           This mode is developed by Drs. Yang Liao, Gordon K Smyth and Wei Shi.
           http://bioinf.wehi.edu.au/featureCounts.


HTseq Modes:

-z 1   HTSeq-Union: Quantify by overlapping.
           If the read (or read pair) overlaps multiple genes, it will be set
           ambiguous. Only reads that overlap one gene will be assigned.

-z 2   HTSeq-Intersection_strict: Quantify by overlapping and intersection.
           This mode requires the assigned gene to cover every base of the read.
           If more than one such genes exist,  the read is set ambiguous.

-z 3   HTSeq-Intersection_nonempty: Quantify by overlapping and intersection.
           This mode does NOT require the assigned gene to cover every base of
           the read, but the gene must cover all sections of the read that overlap
           genes. If more than one such genes exist, the read is set ambiguous.

           Please refer to the HTSeq documentation for more information.
           HTSeq is developed by Dr. Simon Anders at EMBL Heidelberg.
           http://www-huber.embl.de/HTSeq/doc/count.html#count.
           The actual implementation of the HTseq scheme is different in VERSE.


VERSE Modes:

-z 4   Union_strict: A combination of HTSeq-Union and HTSeq-Intersection_strict.
           This mode requires every base of the read to overlap one and only one gene.
           This mode is the most conservative.

-z 5   Cover_length: Quantify by overlapping length comparison.
           After getting a list of overlapping genes, VERSE will calculate the
           overlapping length of each gene and assign the read to the most covered
           gene. One can use --minDifAmbiguous to set the minimum allowed coverage
           difference between the most covered gene and the gene with the second
           highest coverage. If the difference is not large enough the read will
           be set ambiguous.

annotation:
    ## ORFs
    # https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    translation_table: 11
    # when counting reads aligning to ORFs, require at least this many bp
    # overlapping the ORF
    minimum_overlap: 20
    verse_mode: 5

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
        enzyme:
            namemap: /pic/projects/mint/atlas_databases/functional/enzyme/enzyme.db
            fasta: /pic/projects/mint/atlas_databases/functional/enzyme/enzyme.fasta
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
