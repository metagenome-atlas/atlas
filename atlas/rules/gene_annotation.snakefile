def gff_to_gtf(gff_in, gtf_out):
    # orf_re = re.compile(r"ID=(.*?)\;")
    with open(gtf_out, "w") as fh, open(gff_in) as gff:
        for line in gff:
            if line.startswith("##FASTA"): break
            if line.startswith("#"): continue
            # convert:
            # ID=POMFPAEF_00802;inference=ab initio prediction:Prodigal:2.60;
            # to
            # ID POMFPAEF_00802; inference ab initio prediction:Prodigal:2.60;
            toks = line.strip().split("\t")
            toks[-1] = toks[-1].replace("=", " ").replace(";", "; ")
            print(*toks, sep="\t", file=fh)


if config.get("gene_predicter", "prodigal") == "prokka":

    rule run_prokka_annotation:
        input:
            "{sample}/{sample}_contigs.fasta"
        output:
            discrepancy = "{sample}/annotation/prokka/{sample}.err",
            faa = "{sample}/annotation/predicted_genes/{sample}.faa",
            ffn = "{sample}/annotation/predicted_genes/{sample}.ffn",
            fna = "{sample}/annotation/predicted_genes/{sample}.fna",
            fsa = "{sample}/annotation/predicted_genes/{sample}.fsa",
            gff = "{sample}/annotation/predicted_genes/{sample}.gff",
            log = "{sample}/annotation/predicted_genes/{sample}.log",
            tbl = "{sample}/annotation/predicted_genes/{sample}.tbl",
            tsv = "{sample}/annotation/predicted_genes/{sample}.tsv",
            txt = "{sample}/annotation/predicted_genes/{sample}.txt"
        benchmark:
            "logs/benchmarks/prokka/{sample}.txt"
        params:
            outdir = lambda wc, output: os.path.dirname(output.faa),
            kingdom = config.get("prokka_kingdom", PROKKA_KINGDOM)
        conda:
            "%s/prokka.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        shell:
            """prokka --outdir {params.outdir} \
                   --force \
                   --prefix {wildcards.sample} \
                   --locustag {wildcards.sample} \
                   --kingdom {params.kingdom} \
                   --metagenome \
                   --cpus {threads} \
                   {input}
              """


    rule add_contig_metadata:
        input:
            hits = "{sample}/annotation/refseq/{sample}_hits.tsv",
            gff = "{sample}/annotation/predicted_genes/{sample}.gff"
        output:
            temp("{sample}/annotation/refseq/{sample}_hits_plus.tsv")
        shell:
            "atlas munge-blast {input.hits} {input.gff} {output}"


    rule sort_munged_blast_hits:
        # ensure blast hits are grouped by contig, ORF, and then decreasing by bitscore
        input:
            "{sample}/annotation/refseq/{sample}_hits_plus.tsv"
        output:
            "{sample}/annotation/refseq/{sample}_hits_plus_sorted.tsv"
        shell:
            "sort -k1,1 -k2,2 -k13,13rn {input} > {output}"


elif config.get("gene_predicter", "prodigal") == "prodigal":

    rule predict_genes:
        input:
            "{sample}/{sample}_contigs.fasta"
        output:
            fna = "{sample}/annotation/predicted_genes/{sample}.fna",
            faa = "{sample}/annotation/predicted_genes/{sample}.faa",
            gff = "{sample}/annotation/predicted_genes/{sample}.gff"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        log:
            "{sample}/logs/gene_annotation/prodigal.txt"
        benchmark:
            "logs/benchmarks/prokka/{sample}.txt"
        threads:
            1
        shell:
            """
            prodigal -i {input} -o {output.gff} -d {output.fna} \
                -a {output.faa} -p meta -f gff 2> >(tee {log})
            """


    localrules: get_contigs_from_gene_names, rename_genes

    # rule rename_genes:
    #     input:
    #         "{sample}/annotation/predicted_genes/{sample}_ambigous_names.{extension}",
    #         tsv= "{sample}/annotation/predicted_genes/{sample}.tsv"
    #     output:
    #         "{sample}/annotation/predicted_genes/{sample}.{extension}"
    #     wildcard_constraints:
    #         extension= "faa|fna"
    #     run:
    #         with open(output[0],'w') as fout:
    #             with open(input[0]) as fin :
    #                 i=0
    #                 for line in fin:
    #                     if line[0]=='>':
    #                         fout.write(">{sample}_{i}\n".format(i=i,**wildcards))
    #                         i+=1
    #                     else:
    #                         fout.write(line)


    rule get_contigs_from_gene_names:
        input:
            faa = "{sample}/annotation/predicted_genes/{sample}.faa",
        output:
            tsv= "{sample}/annotation/predicted_genes/{sample}.tsv"
        run:
            header = ["gene_id", "Contig", "Gene_nr", "Start", "Stop", "Strand", "Annotation"]
            with open(output.tsv, "w") as tsv:
                tsv.write("\t".join(header) + "\n")
                with open(input.faa) as fin:
                    gene_idx = 0
                    for line in fin:
                        if line[0] == ">":
                            text = line[1:].strip().split(" # ")
                            old_gene_name = text[0]
                            text.remove(old_gene_name)
                            sample, contig_nr, gene_nr = old_gene_name.split("_")
                            tsv.write("{gene_id}\t{sample}_{contig_nr}\t{gene_nr}\t{text}\n".format(
                                text="\t".join(text),
                                gene_id=old_gene_name,
                                i=gene_idx,
                                sample=sample,
                                gene_nr=gene_nr,
                                contig_nr=contig_nr))
                            gene_idx += 1


    rule add_contig_metadata:
        input:
            hits = "{sample}/annotation/refseq/{sample}_hits.tsv",
            tsv = "{sample}/annotation/predicted_genes/{sample}.tsv"
        output:
            "{sample}/annotation/refseq/{sample}_hits_plus_sorted.tsv"
        run:
            import pandas as pd

            tsv = pd.read_table(input.tsv, index_col=0, usecols=[0, 1], squeeze=True)
            hits = pd.read_table(input.hits, header=None)
            contigs = tsv.loc[hits.iloc[:, 0]]
            hits.insert(0, contigs.name, contigs.values)
            hits.sort_values([contigs.name, 0, 2], inplace=True)
            hits.to_csv(output[0], sep="\t", index=False, header=False)


rule convert_gff_to_gtf:
    input:
        "{file}.gff"
    output:
        "{file}.gtf"
    run:
        gff_to_gtf(input[0], output[0])


# TODO: doesn't work for prodigal predictions
rule convert_gff_to_tsv:
    input:
        "{sample}/annotation/predicted_genes/{sample}.gff"
    output:
        "{sample}/annotation/predicted_genes/{sample}_plus.tsv"
    shell:
        """atlas gff2tsv {input} {output}"""


# this rule specifies the more general eggNOG rules
localrules: rename_eggNOG_annotation
rule rename_eggNOG_annotation:
    input:
        "{sample}/annotation/predicted_genes/{sample}.emapper.tsv"
    output:
        "{sample}/annotation/eggNOG.tsv"
    shell:
        "cp {input} {output}"


rule find_counts_per_region:
    input:
        gtf = "{sample}/annotation/predicted_genes/{sample}.gtf",
        bam = "{sample}/sequence_alignment/{sample}.bam"
    output:
        summary = "{sample}/annotation/feature_counts/{sample}_counts.txt.summary",
        counts = "{sample}/annotation/feature_counts/{sample}_counts.txt"
    params:
        min_read_overlap = config.get("minimum_region_overlap", MINIMUM_REGION_OVERLAP),
        paired_only= "-B" if config.get("contig_map_paired_only",CONTIG_MAP_PAIRED_ONLY) else "",
        paired_mode = "-p" if PAIRED_END else "",
        multi_mapping = "-M --fraction" if config.get("contig_count_multi_mapped_reads",CONTIG_COUNT_MULTI_MAPPED_READS) else "--primary",
        feature_counts_allow_overlap = "-O --fraction" if config.get("feature_counts_allow_overlap", FEATURE_COUNTS_ALLOW_OVERLAP) else ""
    log:
        "{sample}/logs/quantify/counts_per_region.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """
        featureCounts \
            --minOverlap {params.min_read_overlap} \
            {params.paired_mode} \
            {params.paired_only} \
            -F GTF \
            -T {threads} \
            {params.multi_mapping} \
            {params.feature_counts_allow_overlap} \
            -t CDS \
            -g ID \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bam} 2> {log}
        """


rule run_diamond_blastp:
    input:
        fasta = "{sample}/annotation/predicted_genes/{sample}.faa",
        db = config["diamond_db"]
    output:
        "{sample}/annotation/refseq/{sample}_hits.tsv"
    benchmark:
        "logs/benchmarks/quantify/run_diamond_blastp/{sample}.txt"
    params:
        tmpdir = "--tmpdir %s" % TMPDIR if TMPDIR else "",
        top_seqs = config.get("diamond_top_seqs", DIAMOND_TOP_SEQS),
        e_value = config.get("diamond_e_value", DIAMOND_E_VALUE),
        min_identity = config.get("diamond_min_identity", DIAMOND_MIN_IDENTITY),
        query_cover = config.get("diamond_query_coverage", DIAMOND_QUERY_COVERAGE),
        gap_open = config.get("diamond_gap_open", DIAMOND_GAP_OPEN),
        gap_extend = config.get("diamond_gap_extend", DIAMOND_GAP_EXTEND),
        block_size = config.get("diamond_block_size", DIAMOND_BLOCK_SIZE),
        index_chunks = config.get("diamond_index_chunks", DIAMOND_INDEX_CHUNKS),
        run_mode = "--more-sensitive" if not config.get("diamond_run_mode", "") == "fast" else ""
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """
        diamond blastp \
            --threads {threads} \
            --outfmt 6 \
            --out {output} \
            --query {input.fasta} \
            --db {input.db} \
            --top {params.top_seqs} \
            --evalue {params.e_value} \
            --id {params.min_identity} \
            --query-cover {params.query_cover} \
            {params.run_mode} \
            --gapopen {params.gap_open} \
            --gapextend {params.gap_extend} \
            {params.tmpdir} \
            --block-size {params.block_size} \
            --index-chunks {params.index_chunks}
        """


rule parse_blastp:
    # assign a taxonomy to contigs using the consensus of the ORF assignments
    input:
        "{sample}/annotation/refseq/{sample}_hits_plus_sorted.tsv"
    output:
        "{sample}/annotation/refseq/{sample}_tax_assignments.tsv"
    params:
        namemap = config["refseq_namemap"],
        treefile = config["refseq_tree"],
        summary_method = config.get("summary_method", SUMMARY_METHOD),
        aggregation_method = config.get("aggregation_method", AGGREGATION_METHOD),
        majority_threshold = config.get("majority_threshold", MAJORITY_THRESHOLD),
        min_identity = config.get("diamond_min_identity", DIAMOND_MIN_IDENTITY),
        min_bitscore = config.get("min_bitscore", MIN_BITSCORE),
        min_length = config.get("min_length", MIN_LENGTH),
        max_evalue = config.get("diamond_e_value", DIAMOND_E_VALUE),
        max_hits = config.get("max_hits", MAX_HITS),
        top_fraction = (100 - config.get("diamond_top_seqs", 5)) * 0.01
    shell:
        """
        atlas refseq \
            --summary-method {params.summary_method} \
            --aggregation-method {params.aggregation_method} \
            --majority-threshold {params.majority_threshold} \
            --min-identity {params.min_identity} \
            --min-bitscore {params.min_bitscore} \
            --min-length {params.min_length} \
            --max-evalue {params.max_evalue} \
            --max-hits {params.max_hits} \
            --top-fraction {params.top_fraction} \
            {input} \
            {params.namemap} \
            {params.treefile} \
            {output}
        """


# TODO: make benchmark
#HIGH throughput : split faa in 1Mio faa chunks for next step
rule eggNOG_homology_search:
    input:
        "%s/eggnog.db" % DBDIR,
        faa = "{folder}/{prefix}.faa",
    output:
        temp("{folder}/{prefix}.emapper.seed_orthologs"),
    params:
        data_dir = DBDIR,
        prefix = "{folder}/{prefix}"
    resources:
        mem = config.get("java_mem", JAVA_MEM)
    threads:
        config["threads"]
    conda:
        "%s/eggNOG.yaml" % CONDAENV
    log:
        "{folder}/logs/{prefix}/eggNOG_homology_search_diamond.log"
    shell:
        """
        emapper.py -m diamond --no_annot --no_file_comments \
            --data_dir {params.data_dir} --cpu {threads} -i {input.faa} \
            -o {params.prefix} --override 2> >(tee {log})
        """


#HIGH throughput : concat emapper.seed_orthologs chunks
# run on single machine
rule eggNOG_annotation:
    input:
        "%s/eggnog.db" % DBDIR,
        seed = rules.eggNOG_homology_search.output
    output:
        temp("{folder}/{prefix}.emapper.annotations")
    params:
        data_dir = DBDIR,
        prefix = "{folder}/{prefix}"
    threads:
        config.get("threads", 1)
    resources:
        mem=config["java_mem"]
    conda:
        "%s/eggNOG.yaml" % CONDAENV
    log:
        "{folder}/logs/{prefix}/eggNOG_annotate_hits_table.log"
    shell:
        """
        emapper.py --annotate_hits_table {input.seed} --no_file_comments \
            --override -o {params.prefix} --cpu {threads} --data_dir {params.data_dir} 2> >(tee {log})
        """


rule add_eggNOG_header:
    input:
        "{folder}/{prefix}.emapper.annotations"
    output:
        "{folder}/{prefix}.emapper.tsv"
    run:
        import pandas as pd

        D = pd.read_table(input[0], header=None)
        D.columns = [
            "query_name",
            "seed_eggNOG_ortholog",
            "seed_ortholog_evalue",
            "seed_ortholog_score",
            "predicted_gene_name",
            "GO_terms",
            "KEGG_KO",
            "BiGG_Reactions",
            "Annotation_tax_scope",
            "Matching_OGs",
            "best_OG|evalue|score",
            "categories",
            "eggNOG_HMM_model_annotation",
        ]

        D.to_csv(output[0],sep="\t",index=False)
