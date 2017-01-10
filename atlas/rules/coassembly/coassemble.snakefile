import json
import os
import re
import sys
from itertools import groupby


def gff_to_gtf(gff_in, gtf_out):
    t = re.compile(r'ID=[0-9]+_([0-9]+);')
    with open(gtf_out, "w") as fh, open(gff_in) as gff:
        for line in gff:
            if line.startswith("#"): continue
            toks = line.strip().split("\t")
            orf = t.findall(toks[-1])[0]
            gene_id = toks[0] + "_" + orf
            toks[-1] = 'gene_id "%s"; %s' % (gene_id, toks[-1])
            print(*toks, sep="\t", file=fh)


def read_fasta(fh):
    """Fasta iterator.

    Accepts file handle of .fasta and yields name and sequence.

    Args:
        fh (file): Open file handle of .fasta file

    Yields:
        tuple: name, sequence

    >>> import os
    >>> from itertools import groupby
    >>> f = open("test.fasta", 'w')
    >>> f.write("@seq1\nACTG")
    >>> f.close()
    >>> f = open("test.fastq")
    >>> for name, seq in read_fastq(f):
            assert name == "seq1"
            assert seq == "ACTG"
    >>> f.close()
    >>> os.remove("test.fasta")
    """
    for header, group in groupby(fh, lambda line: line[0] == '>'):
        if header:
            line = next(group)
            name = line[1:].strip()
        else:
            seq = ''.join(line.strip() for line in group)
            yield name, seq


def print_fasta_record(name, seq, out_handle=sys.stdout, wrap=80):
    """Print a fasta record accounting for line wraps in the sequence.

    Args:
        name (str): name or header for fasta entry
        seq (str): sequence
        out_handle (Optional): open file handle in which to write or stdout
        wrap (Optional[int]) : line width of fasta sequence; None is supported for no wrapping

    """
    print('>', name, sep='', file=out_handle)
    if wrap:
        for i in range(0, len(seq), wrap):
            print(seq[i:i + wrap], file=out_handle)
    else:
        print(seq, file=out_handle)


def split_fasta(fasta, chunk_size=250000):
    chunk_size = int(chunk_size)
    fasta = os.path.expanduser(fasta)
    root, ext = os.path.splitext(fasta)

    file_idx = 0
    with open(fasta) as f:
        for i, (name, seq) in enumerate(read_fasta(f)):
            if i % chunk_size == 0:
                if i == 0:
                    ofh = open("{root}_{idx}{ext}".format(root=root, idx=file_idx, ext=ext), "w")
                    print_fasta_record(name, seq, out_handle=ofh)
                else:
                    ofh.close()
                    file_idx += 1
                    ofh = open("{root}_{idx}{ext}".format(root=root, idx=file_idx, ext=ext), "w")
                    print_fasta_record(name, seq, out_handle=ofh)
            else:
                print_fasta_record(name, seq, out_handle=ofh)
    ofh.close()


rule combine_decontaminated_reads:
    input:
        lambda wc: ["{sample}/quality_control/decontamination/{sample}_pe.fastq.gz".format(sample=sample, norm=NORMALIZATION) for sample in config["samples"]["coassemblies"][wc.coassembly]]
    output:
        "coassemblies/{coassembly}/all_decontamination_reads/{coassembly}_pe.fastq.gz"
    shell:
        "cat {input} > {output}"


rule combine_normalized_reads:
    input:
        lambda wc: ["{sample}/quality_control/{norm}/{sample}_pe.fastq.gz".format(sample=sample, norm=NORMALIZATION) for sample in config["samples"]["coassemblies"][wc.coassembly]]
    output:
        "coassemblies/{coassembly}/all_normalized_reads/{coassembly}_pe.fastq.gz"
    shell:
        "cat {input} > {output}"


rule normalize_combined_reads:
    input:
        "coassemblies/{coassembly}/all_normalized_reads/{coassembly}_pe.fastq.gz"
    output:
        "coassemblies/{coassembly}/quality_control/%s/{coassembly}_pe.fastq.gz" % NORMALIZATION
    params:
        k = config["preprocessing"]["normalization"].get("k", 21),
        t = config["preprocessing"]["normalization"].get("t", 100),
        minkmers = config["preprocessing"]["normalization"].get("minkmers", 15)
    log:
        "coassemblies/{coassembly}/logs/{coassembly}_%s.log" % NORMALIZATION
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbnorm.sh in={input} out={output} k={params.k} t={params.t} \
               minkmers={params.minkmers} prefilter=t threads={threads} 2> {log}"""


if config.get("assembler", "megahit") == "megahit":
    rule coassembly_megahit:
        input:
            "coassemblies/{coassembly}/quality_control/%s/{coassembly}_pe.fastq.gz" % NORMALIZATION
        output:
            temp("coassemblies/{coassembly}/%s/{coassembly}_prefilter.contigs.fa" % ASSEMBLER)
        params:
            memory = config["assembly"].get("memory", 0.90),
            min_count = config["assembly"].get("minimum_count", 2),
            k_min = config["assembly"].get("kmer_min", 21),
            k_max = config["assembly"].get("kmer_max", 121),
            k_step = config["assembly"].get("kmer_step", 20),
            merge_level = config["assembly"].get("merge_level", "20,0.98"),
            prune_level = config["assembly"].get("prune_level", 2),
            low_local_ratio = config["assembly"].get("low_local_ratio", 0.2),
            min_contig_len = config["assembly"].get("minimum_contig_length", 200),
            outdir = lambda wc: "coassemblies/{coassembly}/{assembler}".format(coassembly=wc.coassembly, assembler=ASSEMBLER)
        log:
            "coassemblies/{coassembly}/%s/{coassembly}.log" % ASSEMBLER
        threads:
            config.get("threads", 1)
        shell:
            """{SHPFXM} megahit --num-cpu-threads {threads} --12 {input} --continue \
                   --k-min {params.k_min} --k-max {params.k_max} --k-step {params.k_step} \
                   --out-dir {params.outdir} --out-prefix {wildcards.coassembly}_prefilter \
                   --min-contig-len {params.min_contig_len} --min-count {params.min_count} \
                   --merge-level {params.merge_level} --prune-level {params.prune_level} \
                   --low-local-ratio {params.low_local_ratio}"""


    rule coassembly_rename_megahit_output:
        input:
            "coassemblies/{coassembly}/%s/{coassembly}_prefilter.contigs.fa" % ASSEMBLER
        output:
            "coassemblies/{coassembly}/%s/{coassembly}_prefilter_contigs.fasta" % ASSEMBLER
        shell:
            "{SHPFXS} cp {input} {output}"

else:
    rule coassembly_spades:
        input:
            "coassemblies/{coassembly}/quality_control/%s/{coassembly}_pe.fastq.gz" % NORMALIZATION
        output:
            temp("coassemblies/{coassembly}/%s/contigs.fasta" % ASSEMBLER)
        params:
            # memory = config["assembly"].get("memory", 0.90)
            k = config["assembly"].get("spades_k", "auto"),
            outdir = lambda wc: "coassemblies/{coassembly}/{assembler}".format(coassembly=wc.coassembly, assembler=ASSEMBLER)
        log:
            "coassemblies/{coassembly}/%s/spades.log" % ASSEMBLER
        threads:
            config.get("threads", 1)
        shell:
            """{SHPFXM} spades.py -t {threads} -o {params.outdir} --meta --12 {input}"""


    rule coassembly_rename_spades_output:
        input:
            "coassemblies/{coassembly}/%s/contigs.fasta" % ASSEMBLER
        output:
            "coassemblies/{coassembly}/%s/{coassembly}_prefilter_contigs.fasta" % ASSEMBLER
        shell:
            "{SHPFXS} cp {input} {output}"


rule coassembly_prefilter_stats:
    input:
        "coassemblies/{coassembly}/%s/{coassembly}_prefilter_contigs.fasta" % ASSEMBLER
    output:
        "{coassembly}/%s/stats/prefilter_contig_stats.txt" % ASSEMBLER
    threads:
        1
    shell:
        "{SHPFXS} stats.sh in={input} format=3 > {output}"


rule coassembly_prefilter_contig_coverage:
    input:
        fasta = "{coassembly}/%s/{coassembly}_prefilter_contigs.fasta" % ASSEMBLER,
        fastq = "{coassembly}/quality_control/decontamination/{coassembly}_pe.fastq.gz"
    output:
        bhist = "{coassembly}/%s/stats/prefilter_base_composition.txt" % ASSEMBLER,
        bqhist = "{coassembly}/%s/stats/prefilter_box_quality.txt" % ASSEMBLER,
        mhist = "{coassembly}/%s/stats/prefilter_mutation_rates.txt" % ASSEMBLER,
        statsfile = "{coassembly}/%s/stats/prefilter_mapping_stats.txt" % ASSEMBLER,
        covstats = "{coassembly}/%s/stats/prefilter_coverage_stats.txt" % ASSEMBLER
    log:
        "{coassembly}/logs/dirty_contig_coverage_stats.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbmap.sh nodisk=t ref={input.fasta} in={input.fastq} fast=t \
               threads={threads} bhist={output.bhist} bqhist={output.bqhist} mhist={output.mhist} \
               statsfile={output.statsfile} covstats={output.covstats} 2> {log}"""


rule coassembly_contig_filter:
    input:
        fasta = "{coassembly}/%s/{coassembly}_prefilter_contigs.fasta" % ASSEMBLER,
        covstats = "{coassembly}/%s/stats/prefilter_coverage_stats.txt" % ASSEMBLER
    output:
        fasta = "{coassembly}/%s/{coassembly}_contigs.fasta" % ASSEMBLER,
        removed_names = "{coassembly}/%s/{coassembly}_discarded_contigs.txt" % ASSEMBLER
    params:
        minc = config["assembly"].get("minc", 5),
        minp = config["assembly"].get("minp", 40),
        minr = config["assembly"].get("minr", 0),
        minl = config["assembly"].get("minl", 1),
        trim = config["assembly"].get("trim", 0)
    log:
        "{coassembly}/logs/filter_by_coverage.log"
    threads:
        1
    shell:
        """{SHPFXS} filterbycoverage.sh in={input.fasta} cov={input.covstats} out={output.fasta} \
               outd={output.removed_names} minc={params.minc} minp={params.minp} \
               minr={params.minr} minl={params.minl} trim={params.trim} 2> {log}"""


rule coassembly_postfilter_stats:
   input:
       "{coassembly}/%s/{coassembly}_contigs.fasta" % ASSEMBLER
   output:
       "{coassembly}/%s/stats/final_contig_stats.txt" % ASSEMBLER
   threads:
       1
   shell:
       "{SHPFXS} stats.sh in={input} format=3 > {output}"


rule coassembly_sample_mapping:
    input:
        fasta = "{coassembly}/%s/{coassembly}_contigs.fasta" % ASSEMBLER,
        # this may need to be a function
        fastq = "{sample}/quality_control/decontamination/{sample}_pe.fastq.gz"
    output:
        sam = temp("coassemblies/{coassembly}/annotation/{sample}.sam"),
    log:
        "coassemblies/{coassembly}/logs/{sample}_mapping_stats.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbmap.sh nodisk=t ref={input.fasta} in={input.fastq} trimreaddescriptions=t \
               out={output.sam} mappedonly=t threads={threads} mdtag=t xstag=fs nmtag=t sam=1.3 \
               2> {log}"""


rule coassembly_sam_to_bam:
   input:
       "coassemblies/{coassembly}/annotation/{sample}.sam"
   output:
       "coassemblies/{coassembly}/annotation/{sample}.bam"
   threads:
       config.get("threads", 1)
   shell:
       """{SHPFXM} samtools view -@ {threads} -bSh1 {input} | samtools sort -@ {threads} -T {TMPDIR}/{wildcards.sample}_tmp -o {output} -O bam -"""


rule coassembly_index_bam:
   input:
       "{sample}/annotation/{sample}.bam"
   output:
       "{sample}/annotation/{sample}.bam.bai"
   threads:
       1
   shell:
       "{SHPFXS} samtools index {input}"


rule coassembly_prodigal:
    input:
        "{sample}/%s/{sample}_contigs.fasta" % ASSEMBLER
    output:
        prot = "{sample}/annotation/orfs/{sample}.faa",
        nuc = "{sample}/annotation/orfs/{sample}.fna",
        gff = "{sample}/annotation/orfs/{sample}.gff"
    params:
        g = config["annotation"].get("translation_table", "11")
    threads:
        1
    shell:
        """{SHPFXS} prodigal -i {input} -o {output.gff} -f gff -a {output.prot} -d {output.nuc} \
               -g {params.g} -p meta"""


rule coassembly_gff_to_gtf:
    input:
        "{sample}/annotation/orfs/{sample}.gff"
    output:
        "{sample}/annotation/orfs/{sample}.gtf"
    run:
        gff_to_gtf(input[0], output[0])


rule coassembly_counts_per_region:
    input:
        gtf = "{sample}/annotation/orfs/{sample}.gtf",
        bam = "{sample}/annotation/{sample}.bam",
        bai = "{sample}/annotation/{sample}.bam.bai"
    output:
        summary = "{sample}/annotation/orfs/{sample}.CDS.summary.txt",
        counts = "{sample}/annotation/orfs/{sample}.CDS.txt"
    params:
        min_read_overlap = config["annotation"].get("minimum_overlap", 20)
    log:
        "{sample}/logs/counts_per_region.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} verse -T {threads} --minReadOverlap {params.min_read_overlap} \
               --singleEnd -t CDS -z 1 -a {input.gtf} \
               -o {wildcards.sample}/annotation/orfs/{wildcards.sample} \
               {input.bam} > {log}"""


rule coassembly_split_orfs:
    input:
        faa = "{sample}/annotation/orfs/{sample}.faa"
    output:
        temp(dynamic("{sample}/annotation/orfs/{sample}_{n}.faa"))
    params:
        chunk_size = config["annotation"].get("chunk_size", 250000)
    run:
        split_fasta(input.faa, chunk_size=params.chunk_size)


rule coassembly_diamond_alignments:
    input:
        fasta = "{sample}/annotation/orfs/{sample}_{n}.faa",
        db = lambda wc: config["annotation"]["references"][wc.reference]["dmnd"]
    output:
        temp("{sample}/annotation/{reference}/{sample}_intermediate_{n}.aln")
    params:
        tmpdir = "--tmpdir %s" % TMPDIR if TMPDIR else "",
        top_seqs = lambda wc: config["annotation"]["references"][wc.reference].get("top_seqs", "5"),
        e_value = lambda wc: config["annotation"]["references"][wc.reference].get("e_value", "0.000001"),
        min_identity = lambda wc: config["annotation"]["references"][wc.reference].get("min_identity", "50"),
        query_cover = lambda wc: config["annotation"]["references"][wc.reference].get("query_coverage", "60"),
        gap_open = lambda wc: config["annotation"]["references"][wc.reference].get("gap_open", "11"),
        gap_extend = lambda wc: config["annotation"]["references"][wc.reference].get("gap_extend", "1"),
        block_size = lambda wc: config["annotation"]["references"][wc.reference].get("block_size", "2"),
        index_chunks = lambda wc: config["annotation"]["references"][wc.reference].get("index_chunks", "4"),
        run_mode = lambda wc: "" if config["annotation"]["references"][wc.reference].get("run_mode", "fast") == "fast" else "--more-sensitive"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} diamond blastp --threads {threads} --outfmt 6 --out {output} \
               --query {input.fasta} --db {input.db} --top {params.top_seqs} \
               --evalue {params.e_value} --id {params.min_identity} \
               --query-cover {params.query_cover} {params.run_mode} --gapopen {params.gap_open} \
               --gapextend {params.gap_extend} {params.tmpdir} --block-size {params.block_size} \
               --index-chunks {params.index_chunks}"""


rule merge_alignments:
    input:
        dynamic("{sample}/annotation/{reference}/{sample}_intermediate_{n}.aln")
    output:
        "{sample}/annotation/{reference}/{sample}_hits.tsv"
    shell:
        "{SHPFXS} cat {input} | sort -k1,1 -k12,12rn > {output}"


rule parse_blast:
    input:
        "{sample}/annotation/{reference}/{sample}_hits.tsv"
    output:
        "{sample}/annotation/{reference}/{sample}_assignments.tsv"
    params:
        namemap = lambda wc: config["annotation"]["references"][wc.reference]["namemap"],
        treefile = lambda wc: config["annotation"]["references"][wc.reference].get("tree", ""),
        summary_method = lambda wc: config["annotation"]["references"][wc.reference].get("summary_method", "best"),
        aggregation_method = lambda wc: "--aggregation-method %s" % config["annotation"]["references"][wc.reference].get("aggregation_method", "") if "refseq" in wc.reference else "",
        majority_threshold = lambda wc: "--majority-threshold %f" % config["annotation"]["references"][wc.reference].get("majority_threshold", 0.51) if "refseq" in wc.reference else "",
        min_identity = lambda wc: config["annotation"]["references"][wc.reference].get("min_identity", "50"),
        min_bitscore = lambda wc: config["annotation"]["references"][wc.reference].get("min_bitscore", "0"),
        min_length = lambda wc: config["annotation"]["references"][wc.reference].get("min_length", "60"),
        max_evalue = lambda wc: config["annotation"]["references"][wc.reference].get("max_evalue", "0.000001"),
        max_hits = lambda wc: config["annotation"]["references"][wc.reference].get("max_hits", "10"),
        top_fraction = lambda wc: config["annotation"]["references"][wc.reference].get("top_fraction", "0.50")
    shell:
        """{SHPFXS} atlas {wildcards.reference} \
               --summary-method {params.summary_method} {params.aggregation_method} \
               {params.majority_threshold} --min-identity {params.min_identity} \
               --min-bitscore {params.min_bitscore} --min-length {params.min_length} \
               --max-evalue {params.max_evalue} --max-hits {params.max_hits} \
               --top-fraction {params.top_fraction} {input} {params.namemap} {params.treefile} \
               {output}"""


rule merge_blast:
    input:
        ["{sample}/annotation/%s/{sample}_assignments.tsv" % i for i in list(config["annotation"]["references"].keys())]
    output:
        "{sample}/annotation/{sample}_merged_assignments.tsv"
    shell:
        "{SHPFXS} atlas merge-tables {input} {output}"


rule aggregate_counts:
    input:
        merged = "{sample}/annotation/{sample}_merged_assignments.tsv",
        counts = "{sample}/annotation/orfs/{sample}.CDS.txt"
    output:
        ["{sample}/count_tables/{sample}_%s.tsv" % i for i in TABLES]
    params:
        prefix = lambda wc: "{sample}/count_tables/{sample}".format(sample=wc.sample),
        combos = json.dumps(config["summary_counts"])
    shell:
        """{SHPFXS} atlas counts {params.prefix} {input.merged} \
               {input.counts} '{params.combos}'"""


# compile these counts into a single table with metric, sample count
