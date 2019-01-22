# def gff_to_gtf(gff_in, gtf_out):
#     # orf_re = re.compile(r"ID=(.*?)\;")
#     with open(gtf_out, "w") as fh, open(gff_in) as gff:
#         for line in gff:
#             if line.startswith("##FASTA"): break
#             if line.startswith("#"): continue
#             # convert:
#             # ID=POMFPAEF_00802;inference=ab initio prediction:Prodigal:2.60;
#             # to
#             # ID POMFPAEF_00802; inference ab initio prediction:Prodigal:2.60;
#             toks = line.strip().split("\t")
#             toks[-1] = toks[-1].replace("=", " ").replace(";", "; ")
#             print(*toks, sep="\t", file=fh)
#
# #
# if config.get("run_prokka", False):
#
#     rule run_prokka_annotation:
#         input:
#             "{sample}/{sample}_contigs.fasta"
#         output:
#             discrepancy = "{sample}/annotation/prokka/{sample}.err",
#             faa = "{sample}/annotation/predicted_genes/{sample}.faa",
#             ffn = "{sample}/annotation/predicted_genes/{sample}.ffn",
#             fna = "{sample}/annotation/predicted_genes/{sample}.fna",
#             fsa = "{sample}/annotation/predicted_genes/{sample}.fsa",
#             gff = "{sample}/annotation/predicted_genes/{sample}.gff",
#             log = "{sample}/annotation/predicted_genes/{sample}.log",
#             tbl = "{sample}/annotation/predicted_genes/{sample}.tbl",
#             tsv = "{sample}/annotation/predicted_genes/{sample}_prokka.tsv",
#             txt = "{sample}/annotation/predicted_genes/{sample}.txt"
#         benchmark:
#             "logs/benchmarks/prokka/{sample}.txt"
#         params:
#             outdir = lambda wc, output: os.path.dirname(output.faa),
#             kingdom = config.get("prokka_kingdom", PROKKA_KINGDOM)
#         conda:
#             "%s/prokka.yaml" % CONDAENV
#         threads:
#             config.get("threads", 1)
#         shell:
#             """prokka --outdir {params.outdir} \
#                    --force \
#                    --prefix {wildcards.sample} \
#                    --locustag {wildcards.sample} \
#                    --kingdom {params.kingdom} \
#                    --metagenome \
#                    --cpus {threads} \
#                    {input}
#               """
#
#     # TODO: doesn't work for prodigal predictions
#     rule convert_gff_to_tsv:
#         input:
#             "{sample}/annotation/predicted_genes/{sample}.gff"
#         output:
#             "{sample}/annotation/predicted_genes/{sample}.tsv"
#         shell:
#             """atlas gff2tsv {input} {output}"""
#
#
#     rule add_contig_metadata:
#         input:
#             hits = "{sample}/annotation/refseq/{sample}_hits.tsv",
#             gff = "{sample}/annotation/predicted_genes/{sample}.gff"
#         output:
#             temp("{sample}/annotation/refseq/{sample}_hits_plus.tsv")
#         shell:
#             "atlas munge-blast {input.hits} {input.gff} {output}"
#
#
#
#
# else:





    # localrules: rename_genes
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

##     #
#     # rule add_contig_metadata:
#     #     input:
#     #         hits = "{sample}/annotation/refseq/{sample}_hits.tsv",
#     #         tsv = "{sample}/annotation/predicted_genes/{sample}.tsv"
#     #     output:
#     #         "{sample}/annotation/refseq/{sample}_hits_plus_sorted.tsv"
#     #     run:
#     #         import pandas as pd
#     #
#     #         tsv = pd.read_table(input.tsv, index_col=0, usecols=[0, 1], squeeze=True)
#     #         hits = pd.read_table(input.hits, header=None)
#     #         contigs = tsv.loc[hits.iloc[:, 0]]
#     #         hits.insert(0, contigs.name, contigs.values)
#     #         hits.sort_values([contigs.name, 0, 2], inplace=True)
#     #         hits.to_csv(output[0], sep="\t", index=False, header=False)
#     #
#
# rule convert_gff_to_gtf:
#     input:
#         "{file}.gff"
#     output:
#         "{file}.gtf"
#     run:
#         gff_to_gtf(input[0], output[0])
#
#
#
# no longer used

# rule find_counts_per_region:
#     input:
#         gtf = "{sample}/annotation/predicted_genes/{sample}.gtf",
#         bam = "{sample}/sequence_alignment/{sample}.bam"
#     output:
#         summary = "{sample}/annotation/feature_counts/{sample}_counts.txt.summary",
#         counts = "{sample}/annotation/feature_counts/{sample}_counts.txt"
#     params:
#         min_read_overlap = config.get("minimum_region_overlap", MINIMUM_REGION_OVERLAP),
#         paired_only= "-B" if config.get("contig_map_paired_only",CONTIG_MAP_PAIRED_ONLY) else "",
#         paired_mode = "-p" if PAIRED_END else "",
#         multi_mapping = "-M --fraction" if config.get("contig_count_multi_mapped_reads",CONTIG_COUNT_MULTI_MAPPED_READS) else "--primary",
#         feature_counts_allow_overlap = "-O --fraction" if config.get("feature_counts_allow_overlap", FEATURE_COUNTS_ALLOW_OVERLAP) else ""
#     log:
#         "{sample}/logs/quantify/counts_per_region.log"
#     conda:
#         "%s/required_packages.yaml" % CONDAENV
#     threads:
#         config.get("threads", 1)
#     shell:
#         """
#         featureCounts \
#             --minOverlap {params.min_read_overlap} \
#             {params.paired_mode} \
#             {params.paired_only} \
#             -F GTF \
#             -T {threads} \
#             {params.multi_mapping} \
#             {params.feature_counts_allow_overlap} \
#             -t CDS \
#             -g ID \
#             -a {input.gtf} \
#             -o {output.counts} \
#             {input.bam} 2> {log}
#         """

#
# #### Taxonomy ####
#
#
# rule run_diamond_blastp:
#     input:
#         fasta = "{sample}/annotation/predicted_genes/{sample}.faa",
#         db = config["diamond_db"]
#     output:
#         "{sample}/annotation/refseq/{sample}_hits.tsv"
#     benchmark:
#         "logs/benchmarks/quantify/run_diamond_blastp/{sample}.txt"
#     params:
#         tmpdir = "--tmpdir %s" % TMPDIR if TMPDIR else "",
#         top_seqs = config.get("diamond_top_seqs", DIAMOND_TOP_SEQS),
#         e_value = config.get("diamond_e_value", DIAMOND_E_VALUE),
#         min_identity = config.get("diamond_min_identity", DIAMOND_MIN_IDENTITY),
#         query_cover = config.get("diamond_query_coverage", DIAMOND_QUERY_COVERAGE),
#         gap_open = config.get("diamond_gap_open", DIAMOND_GAP_OPEN),
#         gap_extend = config.get("diamond_gap_extend", DIAMOND_GAP_EXTEND),
#         block_size = config.get("diamond_block_size", DIAMOND_BLOCK_SIZE),
#         index_chunks = config.get("diamond_index_chunks", DIAMOND_INDEX_CHUNKS),
#         run_mode = "--more-sensitive" if not config.get("diamond_run_mode", "") == "fast" else ""
#     conda:
#         "%s/required_packages.yaml" % CONDAENV
#     threads:
#         config.get("threads", 1)
#     shell:
#         """
#         diamond blastp \
#             --threads {threads} \
#             --outfmt 6 \
#             --out {output} \
#             --query {input.fasta} \
#             --db {input.db} \
#             --top {params.top_seqs} \
#             --evalue {params.e_value} \
#             --id {params.min_identity} \
#             --query-cover {params.query_cover} \
#             {params.run_mode} \
#             --gapopen {params.gap_open} \
#             --gapextend {params.gap_extend} \
#             {params.tmpdir} \
#             --block-size {params.block_size} \
#             --index-chunks {params.index_chunks}
#         """
# rule sort_munged_blast_hits:
#     # ensure blast hits are grouped by contig, ORF, and then decreasing by bitscore
#     input:
#         "{sample}/annotation/refseq/{sample}_hits_plus.tsv"
#     output:
#         "{sample}/annotation/refseq/{sample}_hits_plus_sorted.tsv"
#     shell:
#         "sort -k1,1 -k2,2 -k13,13rn {input} > {output}"
#
#
# rule parse_blastp:
#     # assign a taxonomy to contigs using the consensus of the ORF assignments
#     input:
#         "{sample}/annotation/refseq/{sample}_hits_plus_sorted.tsv"
#     output:
#         "{sample}/annotation/refseq/{sample}_tax_assignments.tsv"
#     params:
#         namemap = config["refseq_namemap"],
#         treefile = config["refseq_tree"],
#         summary_method = config.get("summary_method", SUMMARY_METHOD),
#         aggregation_method = config.get("aggregation_method", AGGREGATION_METHOD),
#         majority_threshold = config.get("majority_threshold", MAJORITY_THRESHOLD),
#         min_identity = config.get("diamond_min_identity", DIAMOND_MIN_IDENTITY),
#         min_bitscore = config.get("min_bitscore", MIN_BITSCORE),
#         min_length = config.get("min_length", MIN_LENGTH),
#         max_evalue = config.get("diamond_e_value", DIAMOND_E_VALUE),
#         max_hits = config.get("max_hits", MAX_HITS),
#         top_fraction = (100 - config.get("diamond_top_seqs", 5)) * 0.01
#     shell:
#         """
#         atlas refseq \
#             --summary-method {params.summary_method} \
#             --aggregation-method {params.aggregation_method} \
#             --majority-threshold {params.majority_threshold} \
#             --min-identity {params.min_identity} \
#             --min-bitscore {params.min_bitscore} \
#             --min-length {params.min_length} \
#             --max-evalue {params.max_evalue} \
#             --max-hits {params.max_hits} \
#             --top-fraction {params.top_fraction} \
#             {input} \
#             {params.namemap} \
#             {params.treefile} \
#             {output}
#         """
