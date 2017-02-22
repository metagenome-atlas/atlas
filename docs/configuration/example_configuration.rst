.. _example-configuration:

Example Configuration
=====================

Lines starting with '#' are treated as comments::

    samples:
        # sample ID
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
        coassemblies:
            New-Name:
                - Sample-1
                - Sample-2

    tmpdir: /scratch

    threads: 24

    preprocessing:
        adapters: /databases/adapters.fa

        contamination:
            references:
                rRNA: /databases/silva_rfam_all_rRNAs.fa
                phiX: /databases/phiX174_virus.fa
            ambiguous: best

        normalization:
            k: 19
            t: 100

    assembly:
        kmer_min: 21
        kmer_max: 121
        kmer_step: 20
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
        # count primary read alignments only when getting counts across ORFs
        primary_only: false
        # allow counting of multi-mapped reads when getting counts across ORFs
        multi_mapping: true

        references:
            eggnog:
                namemap: /databases/eggnog.db
                dmnd: /databases/eggnog.dmnd
                block_size: 4
            refseq:
                namemap: /databases/refseq.db
                tree: /databases/refseq.tree
                dmnd: /databases/refseq.dmnd
                block_size: 6
                index_chunks: 1
            enzyme:
                namemap: /databases/enzyme.db
                dmnd: /databases/enzyme.dmnd
                chunk_size: 500000
                top_seqs: 2
                index_chunks: 1
                summary_method: majority
            cazy:
                namemap: /databases/cazy.db
                dmnd: /databases/cazy.dmnd
                chunk_size: 500000
                top_seqs: 2
                index_chunks: 1
                summary_method: majority
            cog:
                namemap: /databases/cog.db
                dmnd: /databases/cog.dmnd
                chunk_size: 500000
                top_seqs: 2
                index_chunks: 1
                summary_method: majority

    summary_counts:
        taxonomy:
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
                - cog_functional_class
                - cog_annotation
            CAZy_EC:
                - cazy_ec
            CAZy_family:
                - cazy_family
            ENZYME:
                - enzyme_name
                - enzyme_ec
            RefSeq:
                - refseq_product
        KO:
            - ko_id
            - ko_gene_symbol
            - ko_product
            - ko_ec
        CAZY_EC:
            - cazy_ec
        COG:
            - cog_id
            - cog_functional_class
            - cog_annotation
        ENZYME:
            - enzyme_name
            - enzyme_ec
        RefSeq:
            - refseq_product
