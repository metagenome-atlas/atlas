import json
import os
import re
import sys
from itertools import groupby


def get_handle_ribosomal_rna_input(wildcards):
    inputs = []
    data_type = config["samples"][wildcards.sample].get("type", "metagenome").lower()
    if data_type == "metagenome":
        return ["{sample}/quality_control/decontamination/{sample}_clean.fastq.gz".format(sample=wildcards.sample),
                "{sample}/quality_control/decontamination/{sample}_rRNA.fastq.gz".format(sample=wildcards.sample)]
    else:
        return ["{sample}/quality_control/decontamination/{sample}_clean.fastq.gz".format(sample=wildcards.sample)]


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


rule quality_filter_reads:
    input:
        lambda wc: config["samples"][wc.sample]["path"]
    output:
        pe = "{sample,^(?!.*coassemblies).*$}/quality_control/quality_filter/{sample}_pe.fastq.gz",
        se = "{sample}/quality_control/quality_filter/{sample}_se.fastq.gz",
        stats = "{sample}/logs/{sample}_quality_filtering_stats.txt"
    params:
        lref = config["preprocessing"]["adapters"],
        rref = config["preprocessing"]["adapters"],
        mink = config["preprocessing"].get("mink", 8),
        trimq = config["preprocessing"].get("minimum_base_quality", 10),
        hdist = config["preprocessing"].get("allowable_kmer_mismatches", 1),
        k = config["preprocessing"].get("reference_kmer_match_length", 27),
        qtrim = config["preprocessing"].get("qtrim", "rl"),
        minlength = config["preprocessing"].get("minimum_passing_read_length", 51),
        minbasefrequency = config["preprocessing"].get("min_base_frequency", 0.05),
        inputs = lambda wc: "in=%s" % config["samples"][wc.sample]["path"][0] if len(config["samples"][wc.sample]["path"]) == 1 else "in=%s in2=%s" % (config["samples"][wc.sample]["path"][0], config["samples"][wc.sample]["path"][1])
    log:
        "{sample}/logs/{sample}_quality_filter_first_pass.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbduk2.sh {params.inputs} out={output.pe} outs={output.se} \
               rref={params.rref} lref={params.lref} mink={params.mink} qout=33 \
               stats={output.stats} hdist={params.hdist} k={params.k} \
               trimq={params.trimq} qtrim={params.qtrim} threads={threads} \
               minlength={params.minlength} minbasefrequency={params.minbasefrequency} \
               overwrite=true 2> {log}"""


rule error_correction:
    input:
        "{sample}/quality_control/quality_filter/{sample}_pe.fastq.gz"
    output:
        "{sample}/quality_control/error_correction/{sample}_pe.fastq.gz"
    log:
        "{sample}/logs/{sample}_error_correction.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} tadpole.sh in={input} out={output} mode=correct threads={threads} \
               ecc=t ecco=t 2> {log}"""


rule decontamination:
    input:
        "{sample}/quality_control/error_correction/{sample}_pe.fastq.gz"
    output:
        dbs = ["{sample}/quality_control/decontamination/{sample}_%s.fastq.gz" % db for db in list(config["preprocessing"]["contamination"]["references"].keys())],
        stats = "{sample}/quality_control/decontamination/{sample}_refstats.txt",
        clean = "{sample}/quality_control/decontamination/{sample}_clean.fastq.gz"
    params:
        refs_in = " ".join(["ref_%s=%s" % (n, fa) for n, fa in config["preprocessing"]["contamination"]["references"].items()]),
        refs_out = lambda wc: " ".join(["out_{ref}={sample}/quality_control/decontamination/{sample}_{ref}.fastq.gz".format(ref=n, sample=wc.sample) for n in list(config["preprocessing"]["contamination"]["references"].keys())]),
        maxindel = config["preprocessing"]["contamination"].get("maxindel", 20),
        minratio = config["preprocessing"]["contamination"].get("minratio", 0.65),
        minhits = config["preprocessing"]["contamination"].get("minhits", 1),
        ambiguous = config["preprocessing"]["contamination"].get("ambiguous", "best"),
        k = config["preprocessing"]["contamination"].get("k", 13)
    log:
        "{sample}/logs/{sample}_decontamination.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbsplit.sh {params.refs_in} in={input} outu={output.clean} \
               {params.refs_out} maxindel={params.maxindel} minratio={params.minratio} \
               minhits={params.minhits} ambiguous={params.ambiguous} refstats={output.stats}\
               threads={threads} k={params.k} local=t 2> {log}"""


rule handle_ribosomal_rna:
    input:
        get_handle_ribosomal_rna_input
    output:
        "{sample}/quality_control/decontamination/{sample}_pe.fastq.gz"
    threads:
        1
    shell:
        "{SHPFXS} cat {input} > {output}"


rule normalization:
    input:
        "{sample}/quality_control/decontamination/{sample}_pe.fastq.gz"
    output:
        "{sample}/quality_control/%s/{sample}_pe.fastq.gz" % NORMALIZATION
    params:
        k = config["preprocessing"]["normalization"].get("k", 21),
        t = config["preprocessing"]["normalization"].get("t", 100),
        minkmers = config["preprocessing"]["normalization"].get("minkmers", 15)
    log:
        "{sample}/logs/{sample}_%s.log" % NORMALIZATION
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbnorm.sh in={input} out={output} k={params.k} t={params.t} \
               minkmers={params.minkmers} prefilter=t threads={threads} 2> {log}"""


if config.get("assembler", "megahit") == "megahit":
    rule megahit:
        input:
            "{sample}/quality_control/%s/{sample}_pe.fastq.gz" % NORMALIZATION
        output:
            temp("{sample}/{assembler}/{sample}_prefilter.contigs.fa")
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
            outdir = lambda wc: "{sample}/{assembler}".format(sample=wc.sample, assembler=wc.assembler)
        log:
            "{sample}/{assembler}/{sample}.log"
        threads:
            config.get("threads", 1)
        shell:
            """{SHPFXM} megahit --num-cpu-threads {threads} --12 {input} --continue \
                   --k-min {params.k_min} --k-max {params.k_max} --k-step {params.k_step} \
                   --out-dir {params.outdir} --out-prefix {wildcards.sample}_prefilter \
                   --min-contig-len {params.min_contig_len} --min-count {params.min_count} \
                   --merge-level {params.merge_level} --prune-level {params.prune_level} \
                   --low-local-ratio {params.low_local_ratio}"""


    rule rename_megahit_output:
        input:
            "{sample}/{assembler}/{sample}_prefilter.contigs.fa"
        output:
            "{sample}/{assembler}/{sample}_prefilter_contigs.fasta"
        shell:
            "{SHPFXS} cp {input} {output}"

else:
    rule spades:
        input:
            "{sample}/quality_control/%s/{sample}_pe.fastq.gz" % NORMALIZATION
        output:
            temp("{sample}/{assembler}/contigs.fasta")
        params:
            # memory = config["assembly"].get("memory", 0.90)
            k = config["assembly"].get("spades_k", "auto"),
            outdir = lambda wc: "{sample}/{assembler}".format(sample=wc.sample, assembler=ASSEMBLER)
        log:
            "{sample}/{assembler}/spades.log"
        threads:
            config.get("threads", 1)
        shell:
            """{SHPFXM} spades.py -t {threads} -o {params.outdir} --meta --12 {input}"""


    rule rename_spades_output:
        input:
            "{sample}/{assembler}/contigs.fasta"
        output:
            "{sample}/{assembler}/{sample}_prefilter_contigs.fasta"
        shell:
            "{SHPFXS} cp {input} {output}"


rule dirty_contigs_stats:
    input:
        "{sample}/{assembler}/{sample}_prefilter_contigs.fasta"
    output:
        "{sample}/{assembler}/stats/prefilter_contig_stats.txt"
    threads:
        1
    shell:
        "{SHPFXS} stats.sh in={input} format=3 > {output}"


rule dirty_contig_coverage_stats:
    input:
        fasta = "{sample}/{assembler}/{sample}_prefilter_contigs.fasta",
        fastq = "{sample}/quality_control/decontamination/{sample}_pe.fastq.gz"
    output:
        bhist = "{sample}/{assembler}/stats/prefilter_base_composition.txt",
        bqhist = "{sample}/{assembler}/stats/prefilter_box_quality.txt",
        mhist = "{sample}/{assembler}/stats/prefilter_mutation_rates.txt",
        statsfile = "{sample}/{assembler}/stats/prefilter_mapping_stats.txt",
        covstats = "{sample}/{assembler}/stats/prefilter_coverage_stats.txt"
    log:
        "{sample}/{assembler}/logs/dirty_contig_coverage_stats.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbmap.sh nodisk=t ref={input.fasta} in={input.fastq} fast=t \
               threads={threads} bhist={output.bhist} bqhist={output.bqhist} mhist={output.mhist} \
               statsfile={output.statsfile} covstats={output.covstats} 2> {log}"""


rule filter_by_coverage:
    input:
        fasta = "{sample}/{assembler}/{sample}_prefilter_contigs.fasta",
        covstats = "{sample}/{assembler}/stats/prefilter_coverage_stats.txt"
    output:
        fasta = "{sample}/{assembler}/{sample}_contigs.fasta",
        removed_names = "{sample}/{assembler}/{sample}_discarded_contigs.txt"
    params:
        minc = config["assembly"].get("minc", 5),
        minp = config["assembly"].get("minp", 40),
        minr = config["assembly"].get("minr", 0),
        minl = config["assembly"].get("minl", 1),
        trim = config["assembly"].get("trim", 0)
    log:
        "{sample}/{assembler}/logs/filter_by_coverage.log"
    threads:
        1
    shell:
        """{SHPFXS} filterbycoverage.sh in={input.fasta} cov={input.covstats} out={output.fasta} \
               outd={output.removed_names} minc={params.minc} minp={params.minp} \
               minr={params.minr} minl={params.minl} trim={params.trim} 2> {log}"""


rule contig_coverage_stats:
    input:
        fasta = "{sample}/{assembler}/{sample}_contigs.fasta",
        fastq = "{sample}/quality_control/decontamination/{sample}_pe.fastq.gz"
    output:
        sam = temp("{sample}/{assembler}/annotation/{sample}.sam"),
        bhist = "{sample}/{assembler}/stats/postfilter_base_composition.txt",
        bqhist = "{sample}/{assembler}/stats/postfilter_box_quality.txt",
        mhist = "{sample}/{assembler}/stats/postfilter_mutation_rates.txt",
        gchist = "{sample}/{assembler}/stats/postfilter_gc_rates.txt",
        statsfile = "{sample}/{assembler}/stats/postfilter_mapping_stats.txt",
        covstats = "{sample}/{assembler}/stats/postfilter_coverage_stats.txt"
    log:
        "{sample}/{assembler}/logs/contig_coverage_stats.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbmap.sh nodisk=t ref={input.fasta} in={input.fastq} trimreaddescriptions=t \
               out={wildcards.sample}/{assembler}/annotation/{wildcards.sample}.sam \
               mappedonly=t threads={threads} bhist={output.bhist} bqhist={output.bqhist} \
               mhist={output.mhist} gchist={output.gchist} statsfile={output.statsfile} \
               covstats={output.covstats} mdtag=t xstag=fs nmtag=t sam=1.3 2> {log}"""


rule sam_to_bam:
   input:
       "{sample}/{assembler}/annotation/{sample}.sam"
   output:
       "{sample}/{assembler}/annotation/{sample}.bam"
   threads:
       config.get("threads", 1)
   shell:
       """{SHPFXM} samtools view -@ {threads} -bSh1 {input} | samtools sort -@ {threads} -T {TMPDIR}/{wildcards.sample}_tmp -o {output} -O bam -"""


rule create_bam_index:
   input:
       "{sample}/{assembler}/annotation/{sample}.bam"
   output:
       "{sample}/{assembler}/annotation/{sample}.bam.bai"
   threads:
       1
   shell:
       "{SHPFXS} samtools index {input}"


rule final_contigs_stats:
   input:
       "{sample}/{assembler}/{sample}_contigs.fasta"
   output:
       "{sample}/{assembler}/stats/final_contig_stats.txt"
   threads:
       1
   shell:
       "{SHPFXS} stats.sh in={input} format=3 > {output}"


rule prodigal_orfs:
    input:
        "{sample}/{assembler}/{sample}_contigs.fasta"
    output:
        prot = "{sample}/{assembler}/annotation/orfs/{sample}.faa",
        nuc = "{sample}/{assembler}/annotation/orfs/{sample}.fna",
        gff = "{sample}/{assembler}/annotation/orfs/{sample}.gff"
    params:
        g = config["annotation"].get("translation_table", "11")
    threads:
        1
    shell:
        """{SHPFXS} prodigal -i {input} -o {output.gff} -f gff -a {output.prot} -d {output.nuc} \
               -g {params.g} -p meta"""


rule gff_to_gtf:
    input:
        "{sample}/{assembler}/annotation/orfs/{sample}.gff"
    output:
        "{sample}/{assembler}/annotation/orfs/{sample}.gtf"
    run:
        gff_to_gtf(input[0], output[0])


rule counts_per_region:
    input:
        gtf = "{sample}/{assembler}/annotation/orfs/{sample}.gtf",
        bam = "{sample}/{assembler}/annotation/{sample}.bam",
        bai = "{sample}/{assembler}/annotation/{sample}.bam.bai"
    output:
        summary = "{sample}/{assembler}/annotation/orfs/{sample}.CDS.summary.txt",
        counts = "{sample}/{assembler}/annotation/orfs/{sample}.CDS.txt"
    params:
        min_read_overlap = config["annotation"].get("minimum_overlap", 20)
    log:
        "{sample}/{assembler}/logs/counts_per_region.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} verse -T {threads} --minReadOverlap {params.min_read_overlap} \
               --singleEnd -t CDS -z 1 -a {input.gtf} \
               -o {wildcards.sample}/{assembler}/annotation/orfs/{wildcards.sample} \
               {input.bam} > {log}"""


rule split:
    input:
        faa = "{sample}/{assembler}/annotation/orfs/{sample}.faa"
    output:
        temp(dynamic("{sample}/{assembler}/annotation/orfs/{sample}_{n}.faa"))
    params:
        chunk_size = config["annotation"].get("chunk_size", 250000)
    run:
        split_fasta(input.faa, chunk_size=params.chunk_size)


rule diamond_alignments:
    input:
        fasta = "{sample}/{assembler}/annotation/orfs/{sample}_{n}.faa",
        db = lambda wc: config["annotation"]["references"][wc.reference]["dmnd"]
    output:
        temp("{sample}/{assembler}/annotation/{reference}/{sample}_intermediate_{n}.aln")
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
        dynamic("{sample}/{assembler}/annotation/{reference}/{sample}_intermediate_{n}.aln")
    output:
        "{sample}/{assembler}/annotation/{reference}/{sample}_hits.tsv"
    shell:
        "{SHPFXS} cat {input} | sort -k1,1 -k12,12rn > {output}"


rule parse_blast:
    input:
        "{sample}/{assembler}/annotation/{reference}/{sample}_hits.tsv"
    output:
        "{sample}/{assembler}/annotation/{reference}/{sample}_assignments.tsv"
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
        ["{sample}/{assembler}/annotation/%s/{sample}_assignments.tsv" % i for i in list(config["annotation"]["references"].keys())]
    output:
        "{sample}/{assembler}/annotation/{sample}_merged_assignments.tsv"
    shell:
        "{SHPFXS} atlas merge-tables {input} {output}"


rule aggregate_counts:
    input:
        merged = "{sample}/{assembler}/annotation/{sample}_merged_assignments.tsv",
        counts = "{sample}/{assembler}/annotation/orfs/{sample}.CDS.txt"
    output:
        ["{sample}/{assembler}/count_tables/{sample}_%s.tsv" % i for i in TABLES]
    params:
        prefix = lambda wc: "{sample}/{assembler}/count_tables/{sample}".format(assembler=ASSEMBLER, sample=wc.sample),
        combos = json.dumps(config["summary_counts"])
    shell:
        """{SHPFXS} atlas counts {params.prefix} {input.merged} \
               {input.counts} '{params.combos}'"""
