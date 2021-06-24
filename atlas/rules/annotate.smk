import re


def gff_to_gtf(gff_in, gtf_out):
    # orf_re = re.compile(r"ID=(.*?)\;")
    with open(gtf_out, "w") as fh, open(gff_in) as gff:
        for line in gff:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            # convert:
            # ID=POMFPAEF_00802;inference=ab initio prediction:Prodigal:2.60;
            # to
            # ID POMFPAEF_00802; inference ab initio prediction:Prodigal:2.60;
            toks = line.strip().split("\t")
            toks[-1] = toks[-1].replace("=", " ").replace(";", "; ")
            print(*toks, sep="\t", file=fh)


def bb_cov_stats_to_maxbin(tsv_in, tsv_out):
    with open(tsv_in) as fi, open(tsv_out, "w") as fo:
        # header
        next(fi)
        for line in fi:
            toks = line.strip().split("\t")
            print(toks[0], toks[1], sep="\t", file=fo)


# TODO this should provide an old name to new name map
rule rename_input_contigs:
    input:
        lambda wc: config["samples"][wc.sample]["fasta"],
    output:
        "{sample}/{sample}_renamed_contigs.fasta",
    conda:
        "%s/required_packages.yaml" % CONDAENV
    shell:
        """rename.sh in={input} out={output} ow=t prefix={wildcards.sample}"""


rule calculate_contigs_stats:
    input:
        "{sample}/{sample}_renamed_contigs.fasta",
    output:
        "{sample}/contig_stats.txt",
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads: 1
    shell:
        "stats.sh in={input} format=3 > {output}"


rule run_prokka_annotation:
    input:
        "{sample}/{sample}_renamed_contigs.fasta",
    output:
        discrepancy="{sample}/prokka/{sample}.err",
        faa="{sample}/prokka/{sample}.faa",
        ffn="{sample}/prokka/{sample}.ffn",
        fna="{sample}/prokka/{sample}.fna",
        fsa="{sample}/prokka/{sample}.fsa",
        gff="{sample}/prokka/{sample}.gff",
        log="{sample}/prokka/{sample}.log",
        tbl="{sample}/prokka/{sample}.tbl",
        tsv="{sample}/prokka/{sample}.tsv",
        txt="{sample}/prokka/{sample}.txt",
    benchmark:
        "logs/benchmarks/annotate/prokka/{sample}.txt"
    params:
        outdir="{sample}/prokka",
        kingdom=config.get("prokka_kingdom", "Bacteria"),
        # metagenome: false
        metagenome="--metagenome" if config.get("metagenome", True) else "",
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads: config.get("threads", 1)
    shell:
        """prokka --outdir {params.outdir} --force --prefix {wildcards.sample} \
        --locustag {wildcards.sample} --kingdom {params.kingdom} {params.metagenome} \
        --cpus {threads} {input}"""


rule update_prokka_tsv:
    input:
        "{sample}/prokka/{sample}.gff",
    output:
        "{sample}/prokka/{sample}_fixed.tsv",
    shell:
        """atlas gff2tsv {input} {output}"""


rule convert_gff_to_gtf:
    input:
        "{sample}/prokka/{sample}.gff",
    output:
        "{sample}/prokka/{sample}.gtf",
    run:
        gff_to_gtf(input[0], output[0])


rule run_diamond_blastp:
    input:
        fasta="{sample}/prokka/{sample}.faa",
        db=config["diamond_db"],
    output:
        "{sample}/refseq/{sample}_hits.tsv",
    benchmark:
        "logs/benchmarks/annotate/diamond_alignments/{sample}.txt"
    params:
        tmpdir="--tmpdir %s" % TMPDIR if TMPDIR else "",
        top_seqs=config.get("diamond_top_seqs", 2),
        e_value=config.get("diamond_e_value", 0.000001),
        min_identity=config.get("diamond_min_identity", 50),
        query_cover=config.get("diamond_query_coverage", 60),
        gap_open=config.get("diamond_gap_open", 11),
        gap_extend=config.get("diamond_gap_extend", 1),
        block_size=config.get("diamond_block_size", 2),
        index_chunks=config.get("diamond_index_chunks", 4),
        run_mode=(
            "--more-sensitive"
            if not config.get("diamond_run_mode", "") == "fast"
            else ""
        ),
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads: config.get("threads", 1)
    shell:
        """diamond blastp \
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
        --index-chunks {params.index_chunks}"""


rule add_contig_metadata_to_blast_hits:
    input:
        hits="{sample}/refseq/{sample}_hits.tsv",
        gff="{sample}/prokka/{sample}.gff",
    output:
        temp("{sample}/refseq/{sample}_hits_plus.tsv"),
    shell:
        "atlas munge-blast {input.hits} {input.gff} {output}"


rule sort_munged_blast_hits:
    # ensure blast hits are grouped by contig, ORF, and then decreasing by bitscore
    input:
        "{sample}/refseq/{sample}_hits_plus.tsv",
    output:
        "{sample}/refseq/{sample}_hits_plus_sorted.tsv",
    shell:
        "sort -k1,1 -k2,2 -k13,13rn {input} > {output}"


rule parse_blastp:
    # assign a taxonomy to contigs using the consensus of the ORF assignments
    input:
        "{sample}/refseq/{sample}_hits_plus_sorted.tsv",
    output:
        "{sample}/refseq/{sample}_tax_assignments.tsv",
    params:
        namemap=config["refseq_namemap"],
        treefile=config["refseq_tree"],
        summary_method=config.get("summary_method", "lca"),
        aggregation_method=config.get("aggregation_method", "lca-majority"),
        majority_threshold=config.get("majority_threshold", 0.51),
        min_identity=config.get("diamond_min_identity", 50),
        min_bitscore=config.get("min_bitscore", 0),
        min_length=config.get("min_length", 20),
        max_evalue=config.get("diamond_e_value", 0.000001),
        max_hits=config.get("max_hits", 100),
        top_fraction=(100 - config.get("diamond_top_seqs", 5)) * 0.01,
    shell:
        """atlas refseq --summary-method {params.summary_method} \
        --aggregation-method {params.aggregation_method} \
        --majority-threshold {params.majority_threshold} \
        --min-identity {params.min_identity} \
        --min-bitscore {params.min_bitscore} \
        --min-length {params.min_length} \
        --max-evalue {params.max_evalue} \
        --max-hits {params.max_hits} \
        --top-fraction {params.top_fraction} \
        {input} {params.namemap} {params.treefile} {output}"""


rule align_reads_to_renamed_contigs:
    input:
        fasta="{sample}/{sample}_renamed_contigs.fasta",
        fastq=lambda wc: config["samples"][wc.sample]["fastq"],
    output:
        sam=temp("{sample}/sequence_alignment/{sample}.sam"),
        bhist="{sample}/contig_stats/base_composition.txt",
        bqhist="{sample}/contig_stats/box_quality.txt",
        mhist="{sample}/contig_stats/mutation_rates.txt",
        gchist="{sample}/contig_stats/gc_rates.txt",
        statsfile="{sample}/contig_stats/mapping_stats.txt",
        covstats="{sample}/contig_stats/coverage_stats.txt",
    benchmark:
        "logs/benchmarks/annotate/bbmap_alignment/{sample}.txt"
    params:
        interleaved=(
            lambda wc: "t"
            if config["samples"][wc.sample].get("paired", True)
            and len(config["samples"][wc.sample]["fastq"]) == 1
            else "f"
        ),
        inputs=(
            lambda wc: "in=%s" % config["samples"][wc.sample]["fastq"][0]
            if len(config["samples"][wc.sample]["fastq"]) == 1
            else "in=%s in2=%s"
            % (
                config["samples"][wc.sample]["fastq"][0],
                config["samples"][wc.sample]["fastq"][1],
            )
        ),
        maxsites=config.get("maximum_counted_map_sites", 10),
    log:
        "{sample}/annotate/sequence_alignment/align_reads.log",
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads: config.get("threads", 1)
    shell:
        """bbmap.sh nodisk=t ref={input.fasta} {params.inputs} trimreaddescriptions=t \
        out={output.sam} mappedonly=t threads={threads} bhist={output.bhist} \
        bqhist={output.bqhist} mhist={output.mhist} gchist={output.gchist} \
        statsfile={output.statsfile} covstats={output.covstats} mdtag=t xstag=fs nmtag=t \
        sam=1.3 local=t ambiguous=all interleaved={params.interleaved} secondary=t \
        maxsites={params.maxsites} 2> {log}"""


rule convert_sam_to_bam:
    input:
        "{file}.sam",
    output:
        "{file}.bam",
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads: config.get("threads", 1)
    resources:
        mem=config.get("threads", 1),
    shell:
        """samtools view \
        -m 1G \
        -@ {threads} \
        -bSh1 {input} | samtools sort \
                            -m 1G \
                            -@ {threads} \
                            -T {wildcards.file}_tmp \
                            -o {output} \
                            -O bam -
        """


rule create_bam_index:
    input:
        "{file}.bam",
    output:
        "{file}.bam.bai",
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads: 1
    shell:
        "samtools index {input}"


rule remove_pcr_duplicates:
    input:
        bam="{sample}/sequence_alignment/{sample}.bam",
        bai="{sample}/sequence_alignment/{sample}.bam.bai",
    output:
        bam="{sample}/sequence_alignment/{sample}_markdup.bam",
        txt="{sample}/sequence_alignment/{sample}_markdup_metrics.txt",
    benchmark:
        "logs/benchmarks/annotate/picard_mark_duplicates/{sample}.txt"
    resources:
        mem=int(config.get("java_mem", "32")),
    log:
        "{sample}/annotate/sequence_alignment/remove_pcr_duplicates.log",
    conda:
        "%s/required_packages.yaml" % CONDAENV
    shell:
        """picard MarkDuplicates -Xmx{resources.mem}G INPUT={input.bam} \
        OUTPUT={output.bam} METRICS_FILE={output.txt} ASSUME_SORT_ORDER=coordinate \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=TRUE \
        VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE 2> {log}"""


rule counts_per_region:
    input:
        gtf="{sample}/prokka/{sample}.gtf",
        bam="{sample}/sequence_alignment/{sample}_markdup.bam",
    output:
        summary="{sample}/feature_counts/{sample}_counts.txt.summary",
        counts="{sample}/feature_counts/{sample}_counts.txt",
    params:
        min_read_overlap=config.get("minimum_overlap", 1),
        paired_mode=(
            lambda wc: "-p" if config["samples"][wc.sample].get("paired", True) else ""
        ),
        multi_mapping="-M" if config.get("count_multi_mapped_reads", False) else "",
        primary_only="--primary" if config.get("primary_only", False) else "",
    log:
        "{sample}/annotate/feature_counts/feature_counts.log",
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads: config.get("threads", 1)
    shell:
        """featureCounts {params.paired_mode} -F gtf -T {threads} \
        {params.multi_mapping} -t CDS -g ID -a {input.gtf} -o {output.counts} \
        {input.bam} 2> {log}"""


if config.get("quantification"):

    rule merge_quantification_tables:
        input:
            prokka="{sample}/prokka/{sample}_fixed.tsv",
            refseq="{sample}/refseq/{sample}_tax_assignments.tsv",
            counts="{sample}/feature_counts/{sample}_counts.txt",
        output:
            "{sample}_annotations.txt",
        shell:
            "atlas merge-tables --counts {input.counts} {input.prokka} {input.refseq} {output}"


else:

    rule merge_quantification_tables:
        input:
            prokka="{sample}/prokka/{sample}_fixed.tsv",
            refseq="{sample}/refseq/{sample}_tax_assignments.tsv",
        output:
            "{sample}_annotations.txt",
        shell:
            "atlas merge-tables {input.prokka} {input.refseq} {output}"


rule assembly_report:
    input:
        contig_stats=expand("{sample}/contig_stats.txt", sample=SAMPLES),
        gene_tables=expand("{sample}/prokka/{sample}_fixed.tsv", sample=SAMPLES),
        mapping_log_files=expand(
            "{sample}/sequence_alignment/align_reads.log", sample=SAMPLES
        ),
    output:
        report="reports/assembly_report.html",
        combined_contig_stats="stats/combined_contig_stats.tsv",
    params:
        samples=SAMPLES,
    conda:
        "%s/report.yaml" % CONDAENV
    script:
        "../report/assembly_report.py"
