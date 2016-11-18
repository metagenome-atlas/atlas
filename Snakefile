import os
from glob import glob
from subprocess import check_output


def read_count(fastq):
    total = 0
    count_file = fastq + '.count'
    if os.path.exists(fastq) and os.path.getsize(fastq) > 100:
        if not os.path.exists(count_file):
            check_output("awk '{n++}END{print n/4}' %s > %s" % (fastq, fastq + '.count'), shell=True)
        with open(count_file) as fh:
            for line in fh:
                total = int(line.strip())
                break
    return total


def get_samples(path, coverage_cutoff=1000):
    """Expecting files with .fastq and a naming convention like <sample>_R1.fastq and
    <sample>_R2.fastq.
    """
    samples = set()
    for f in os.listdir(path):
        if f.endswith("fastq") and ("_r1" in f or "_R1" in f):
            if read_count(os.path.join(path, f)) > coverage_cutoff:
                samples.add(f.partition(".")[0].partition("_")[0])
    return samples


def pattern_search(path, patterns):
    """Grab files under path that match pattern."""
    files = []
    for pattern in patterns:
        files.extend(glob("%s/%s" % (path, pattern)))
    return [os.path.basename(i) for i in files]


EID = config['eid']
SAMPLES = get_samples(os.path.join("data", EID), 200)
DECON_DBS = list(config["contamination_filtering"]["references"].keys())


rule all:
    input:
        expand("results/{eid}/decon/{sample}_{decon_dbs}.fastq.gz", eid=EID, sample=SAMPLES, decon_dbs=DECON_DBS),
        expand("results/{eid}/decon/{sample}_refstats.txt", eid=EID, sample=SAMPLES),
        expand("results/{eid}/fastqc/{sample}_final_fastqc.zip", eid=EID, sample=SAMPLES),
        expand("results/{eid}/fastqc/{sample}_final_fastqc.html", eid=EID, sample=SAMPLES),
        expand("results/{eid}/assembly/{sample}/{sample}_length_pass.fa", eid=EID, sample=SAMPLES),
        expand("results/{eid}/annotation/orfs/{sample}_length_pass.faa", eid=EID, sample=SAMPLES),
        expand("results/{eid}/annotation/orfs/{sample}_length_pass.CDS.txt", eid=EID, sample=SAMPLES),
        expand("results/{eid}/annotation/{reference}/{sample}_hits.tsv", eid=EID, reference=config["annotation"]["references"].split(","), sample=SAMPLES)


rule quality_filter_reads:
    input:
        r1 = "data/{eid}/{sample}_R1.fastq",
        r2 = "data/{eid}/{sample}_R2.fastq"
    output:
        r1 = "results/{eid}/quality_filter/{sample}_R1.fastq",
        r2 = "results/{eid}/quality_filter/{sample}_R2.fastq",
        stats = "results/{eid}/logs/{sample}_quality_filtering_stats.txt"
    params:
        lref = config['filtering']['adapters'],
        rref = config['filtering']['adapters'],
        mink = config['filtering']['mink'],
        trimq = config['filtering']['minimum_base_quality'],
        hdist = config['filtering']['allowable_kmer_mismatches'],
        k = config['filtering']['reference_kmer_match_length'],
        qtrim = "rl",
        minlength = config['filtering']['minimum_passing_read_length']
    threads:
        config["threads"]["medium"]
    shell:
        """bbduk2.sh -Xmx8g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} \
               rref={params.rref} lref={params.lref} mink={params.mink} \
               stats={output.stats} hdist={params.hdist} k={params.k} \
               trimq={params.trimq} qtrim={params.qtrim} threads={threads} \
               minlength={params.minlength} overwrite=true"""


rule join_reads:
    input:
        r1 = "results/{eid}/quality_filter/{sample}_R1.fastq",
        r2 = "results/{eid}/quality_filter/{sample}_R2.fastq"
    output:
        joined = "results/{eid}/joined/{sample}.extendedFrags.fastq",
        hist = "results/{eid}/joined/{sample}.hist",
        failed_r1 = "results/{eid}/joined/{sample}.notCombined_1.fastq",
        failed_r2 = "results/{eid}/joined/{sample}.notCombined_2.fastq"
    message:
        "Joining reads using `flash`"
    shadow:
        "shallow"
    params:
        output_dir = lambda wildcards: "results/%s/joined/" % wildcards.eid,
        min_overlap = config['merging']['minimum_overlap'],
        max_overlap = config['merging']['maximum_overlap'],
        max_mismatch_density = config['merging']['maximum_mismatch_density'],
        phred_offset = config['phred_offset']
    log:
        "results/{eid}/logs/{sample}_flash.log"
    threads:
        24
    shell:
        """flash {input.r1} {input.r2} --min-overlap {params.min_overlap} \
               --max-overlap {params.max_overlap} --max-mismatch-density {params.max_mismatch_density} \
               --phred-offset {params.phred_offset} --output-prefix {wildcards.sample} \
               --output-directory {params.output_dir} --threads {threads}"""


rule concatenate_joined_reads:
    input:
        joined = "results/{eid}/joined/{sample}.extendedFrags.fastq",
        failed_r1 = "results/{eid}/joined/{sample}.notCombined_1.fastq",
        failed_r2 = "results/{eid}/joined/{sample}.notCombined_2.fastq"
    output:
        "results/{eid}/joined/{sample}_joined.fastq"
    shell:
        "cat {input.joined} {input.failed_r1} {input.failed_r2} > {output}"


rule error_correction:
    input:
        "results/{eid}/joined/{sample}_joined.fastq"
    output:
        "results/{eid}/joined/{sample}_corrected.fastq.gz"
    threads:
        config["threads"]["large"]
    shell:
        "tadpole.sh in={input} out={output} mode=correct threads={threads}"


if config["qual_method"] == "expected_error":
    rule subset_reads_by_quality:
        input: "results/{eid}/joined/{sample}_corrected.fastq.gz"
        output: "results/{eid}/quality_filter/{sample}_filtered.fastq"
        params:
            phred = config.get("phred_offset", 33),
            maxee = config["filtering"].get("maximum_expected_error", 2),
            maxns = config["filtering"].get("maxns", 3)
        threads:
            1
        shell: """vsearch --fastq_filter {input} --fastqout {output} --fastq_ascii {params.phred} \
                      --fastq_maxee {params.maxee} --fastq_maxns {params.maxns}"""
else:
    rule subset_reads_by_quality:
        input:
            "results/{eid}/joined/{sample}_corrected.fastq.gz"
        output:
            "results/{eid}/quality_filter/{sample}_filtered.fastq"
        params:
            adapter_clip = "" if not config["filtering"].get("adapters", "") else "ILLUMINACLIP:%s:%s" % (config["filtering"]["adapters"], config["filtering"].get("adapter_clip", "2:30:10")),
            window_size_qual = "" if not config["filtering"].get("window_size_quality", "") else "SLIDINGWINDOW:%s" % config["filtering"]["window_size_quality"],
            leading = "" if not config["filtering"].get("leading", 0) else "LEADING:%s" % config["filtering"]["leading"],
            trailing = "" if not config["filtering"].get("trailing", 0) else "TRAILING:%s" % config["filtering"]["trailing"],
            crop = "" if not config["filtering"].get("crop", 0) else "CROP:%s" % config["filtering"]["crop"],
            headcrop = "" if not config["filtering"].get("headcrop", 0) else "HEADCROP:%s" % config["filtering"]["headcrop"],
            minlen = "MINLEN:%s" % config["filtering"]["minimum_passing_read_length"]
        threads:
            config["threads"]["large"]
        shell:
            """trimmomatic SE -threads {threads} {input} {output} {params.adapter_clip} \
                   {params.leading} {params.trailing} {params.window_size_qual} {params.minlen}"""


rule decontaminate_joined:
    input:
        "results/{eid}/quality_filter/{sample}_filtered.fastq"
    output:
        dbs = ["results/{eid}/decon/{sample}_%s.fastq.gz" % db for db in DECON_DBS],
        stats = "results/{eid}/decon/{sample}_refstats.txt",
        clean = "results/{eid}/decon/{sample}_clean.fastq.gz"
    params:
        refs_in = " ".join(["ref_%s=%s" % (n, fa) for n, fa in config['contamination_filtering']['references'].items()]),
        refs_out = lambda wildcards: " ".join(["out_%s=results/%s/decon/%s_%s.fastq.gz" % (n, wildcards.eid, wildcards.sample, n) for n in list(config['contamination_filtering']['references'].keys())]),
        path = "databases/contaminant/",
        maxindel = config['contamination_filtering'].get('maxindel', 20),
        minratio = config['contamination_filtering'].get('minratio', 0.65),
        minhits = config['contamination_filtering'].get('minhits', 1),
        ambiguous = config['contamination_filtering'].get('ambiguous', "best"),
        k = config["contamination_filtering"].get("k", 15)
    threads:
        config["threads"]["large"]
    shell:
        """bbsplit.sh {params.refs_in} path={params.path} in={input} outu={output.clean} \
               {params.refs_out} maxindel={params.maxindel} minratio={params.minratio} \
               minhits={params.minhits} ambiguous={params.ambiguous} refstats={output.stats}\
               threads={threads} k={params.k} local=t"""


if config['data_type'] == "metatranscriptome":
    rule ribosomal_rna:
        input:
            "results/{eid}/decon/{sample}_clean.fastq.gz"
        output:
            "results/{eid}/decon/{sample}_final.fastq.gz"
        shell:
            "cp {input} {output}"
else:
    rule ribosomal_rna:
        input:
            clean = "results/{eid}/decon/{sample}_clean.fastq.gz",
            rrna = "results/{eid}/decon/{sample}_rRNA.fastq.gz"
        output:
            "results/{eid}/decon/{sample}_final.fastq.gz"
        shell:
            "cat {input.clean} {input.rrna} > {output}"


rule fastqc:
    input:
        "results/{eid}/decon/{sample}_final.fastq.gz"
    output:
        "results/{eid}/fastqc/{sample}_final_fastqc.zip",
        "results/{eid}/fastqc/{sample}_final_fastqc.html"
    params:
        output_dir = lambda wildcards: "results/{eid}/fastqc/".format(eid=wildcards.eid)
    threads:
        config["threads"]["large"]
    shell:
        "fastqc -t {threads} -f fastq -o {params.output_dir} {input}"


rule megahit_assembly:
    input:
        "results/{eid}/decon/{sample}_final.fastq.gz"
    output:
        "results/{eid}/assembly/{sample}/{sample}.contigs.fa"
    params:
        memory = config['assembly']['memory'],
        min_count = config['assembly']['minimum_count'],
        k_min = config['assembly']['kmer_min'],
        k_max = config['assembly']['kmer_max'],
        k_step = config['assembly']['kmer_step'],
        merge_level = config['assembly']['merge_level'],
        prune_level = config['assembly']['prune_level'],
        low_local_ratio = config['assembly']['low_local_ratio'],
        min_contig_len = config['assembly']['minimum_contig_length'],
        outdir = lambda wildcards: "results/%s/assembly/%s" % (wildcards.eid, wildcards.sample)
    threads:
        config["threads"]["large"]
    log:
        "results/{eid}/assembly/{sample}/{sample}.log"
    shell:
        """megahit --num-cpu-threads {threads} --read {input} --continue \
               --k-min {params.k_min} --k-max {params.k_max} --k-step {params.k_step} \
               --out-dir {params.outdir} --out-prefix {wildcards.sample} \
               --min-contig-len {params.min_contig_len} --min-count {params.min_count} \
               --merge-level {params.merge_level} --prune-level {params.prune_level} \
               --low-local-ratio {params.low_local_ratio}"""


rule length_filter:
    input:
        "results/{eid}/assembly/{sample}/{sample}.contigs.fa"
    output:
        passing = "results/{eid}/assembly/{sample}/{sample}_length_pass.fa",
        failing = "results/{eid}/assembly/{sample}/{sample}_length_fail.fa"
    params:
        min_contig_length = config['assembly']['filtered_contig_length']
    threads:
        1
    shell:
        """python scripts/fastx.py length-filter --min-length {params.min_contig_length} \
               {input} {output.passing} {output.failing}"""


# rule assembly_stats:


rule build_assembly_index:
    input:
        "results/{eid}/assembly/{sample}/{sample}_length_pass.fa"
    output:
        "results/{eid}/assembly/{sample}/{sample}_length_pass.fa.amb",
        "results/{eid}/assembly/{sample}/{sample}_length_pass.fa.ann",
        "results/{eid}/assembly/{sample}/{sample}_length_pass.fa.bwt",
        "results/{eid}/assembly/{sample}/{sample}_length_pass.fa.pac",
        "results/{eid}/assembly/{sample}/{sample}_length_pass.fa.sa"
    shell:
        "bwa index {input}"


rule align_reads_to_assembly:
    input:
        fastq = "results/{eid}/decon/{sample}_final.fastq.gz",
        ref = "results/{eid}/assembly/{sample}/{sample}_length_pass.fa",
        idx = "results/{eid}/assembly/{sample}/{sample}_length_pass.fa.bwt"
    output:
        "results/{eid}/aligned_reads/{sample}.bam"
    threads:
        config["threads"]["large"]
    shell:
        """bwa mem -t {threads} -L 1,1 {input.ref} {input.fastq} \
               | samtools view -@ {threads} -bS - \
               | samtools sort -@ {threads} -T {wildcards.sample} -o {output} -O bam -"""


rule prodigal_orfs:
    input:
        "results/{eid}/assembly/{sample}/{sample}_length_pass.fa"
    output:
        prot = "results/{eid}/annotation/orfs/{sample}_length_pass.faa",
        nuc = "results/{eid}/annotation/orfs/{sample}_length_pass.fna",
        gff = "results/{eid}/annotation/orfs/{sample}_length_pass.gff"
    params:
        g = config['annotation']['translation_table']
    shell:
        """prodigal -i {input} -o {output.gff} -f gff -a {output.prot} -d {output.nuc} \
               -g {params.g} -p meta"""


rule gff_to_gtf:
    input:
        "results/{eid}/annotation/orfs/{sample}_length_pass.gff"
    output:
        "results/{eid}/annotation/orfs/{sample}_length_pass.gtf"
    run:
        import re
        t = re.compile(r'ID=[0-9]+_([0-9]+);')
        with open(output[0], "w") as fh, open(input[0]) as gff:
            print("##gff-version  3", file=fh)
            for line in gff:
                if line.startswith("#"): continue
                toks = line.strip().split("\t")
                orf = t.findall(toks[-1])[0]
                gene_id = toks[0] + "_" + orf
                toks[-1] = toks[-1] + "gene_id " + gene_id + ";"
                print(*toks, sep="\t", file=fh)


rule counts_per_region:
    input:
        gtf = "results/{eid}/annotation/orfs/{sample}_length_pass.gtf",
        bam = "results/{eid}/aligned_reads/{sample}.bam"
    output:
        summary = "results/{eid}/annotation/orfs/{sample}_length_pass.CDS.summary.txt",
        counts = "results/{eid}/annotation/orfs/{sample}_length_pass.CDS.txt"
    params:
        min_read_overlap = config['annotation']['minimum_overlap']
    threads:
        config["threads"]["medium"]
    shell:
        """verse --multithreadDecompress -T {threads} --minReadOverlap {params.min_read_overlap} \
               --singleEnd -t CDS -z 5 -a {input.gtf} \
               -o results/{wildcards.eid}/annotation/orfs/{wildcards.sample}_length_pass \
               {input.bam} 2> /dev/null"""


rule build_dmnd_database:
    input:
        "databases/{reference}.fasta"
    output:
        "databases/{reference}.fasta.dmnd"
    threads:
        config["threads"]["large"]
    shell:
        "diamond makedb --threads {threads} --in {input} --db {input}"


rule split:
    input:
        "results/{eid}/annotation/orfs/{sample}_length_pass.faa"
    output:
        temp(dynamic("results/{eid}/annotation/orfs/{sample}_length_pass_{n}.faa"))
    params:
        chunk_size = config["annotation"].get("chunk_size", "250000")
    shell:
        "python scripts/fastx.py split-fasta --chunk-size {params.chunk_size} {input}"


rule diamond_alignments:
    input:
        fasta = "results/{eid}/annotation/orfs/{sample}_length_pass_{n}.faa",
        db = "databases/{reference}.fasta.dmnd"
    output:
        temp("results/{eid}/annotation/{reference}/{sample}_intermediate_{n}.aln")
    params:
        tmpdir = "--tmpdir %s" % config.get("temporary_directory", "") if config.get("temporary_directory", "") else "",
        top_seqs = config["annotation"].get("top_seqs", "5"),
        e_value = config["annotation"].get("e_value", "0.000001"),
        min_identity = config["annotation"].get("min_identity", "50"),
        query_cover = config["annotation"].get("query_coverage", "60"),
        gap_open = config["annotation"].get("gap_open", "11"),
        gap_extend = config["annotation"].get("gap_extend", "1"),
        block_size = config["annotation"].get("block_size", "2"),
        index_chunks = config["annotation"].get("index_chunks", "4")
    threads:
        config["threads"]["large"]
    shell:
        """diamond blastp --threads {threads} --outfmt 6 --out {output} \
               --query {input.fasta} --db {input.db} --top {params.max_top_seqs} \
               --evalue {params.e_value} --id {params.min_identity} \
               --query-cover {params.query_cover} --more-sensitive --gapopen {params.gap_open} \
               --gapextend {params.gap_extend} {params.tmpdir} --block-size {params.block_size} \
               --index-chunks {params.index_chunks}"""


rule merge_alignments:
    input:
        dynamic("results/{eid}/annotation/{reference}/{sample}_intermediate_{n}.aln")
    output:
        "results/{eid}/annotation/{reference}/{sample}_hits.tsv"
    shell:
        "cat {input} > {output}"


# rule parse_blast:
#     input:
#         "results/{eid}/annotation/{reference}/{sample}_hits.tsv"
#     output:
#         "results/{eid}/annotation/{reference}/{sample}_assignments.tsv"
#     params:
#
#     shell:
#         """python scripts/blast2assignment.py {params.} """
#
#
# rule merge_blast:
#     input:
#         "results/{eid}/annotation/{reference}/{sample}_assignments.tsv"
#     output:
#         "results/{eid}/annotation/{sample}_merged_assignments.tsv"
#     shell:
#         # input list...
#         # script will need to sniff headers or even if it matters
#         """python scripts/blast2assignment.py merge-tables {input}"""


# rule run_maxbin:
#     input:
#         reads = rules.filter_contaminants.output
#         contigs = rules.assemble.output
#     output:
#         bins = "output/{eid}/binning/{sample}.fasta",
#         abundance = "output/{eid}/binning/{sample}.abund1",
#         log = "output/{eid}/binning/{sample}.log",
#         marker = "output/{eid}/binning/{sample}.marker",
#         summary = "output/{eid}/binning/{sample}.summary",
#         tooshort = "output/{eid}/binning/{sample}.tooshort",
#         noclass = "output/{eid}/binning/{sample}.noclass"
#     params:
#         min_contig_len = config['binning']['minimum_contig_length'],
#         max_iteration = config['binning']['maximum_iterations'],
#         prob_threshold = config['binning']['probability_threshold'],
#         markerset = config['binning']['marker_set']
#     threads:
#         config['binning']['threads']
#     shell:
#         """run_MaxBin.pl -contig {input.contigs} -out {output} -reads {input.reads} \
#         -min_contig_length {params.min_contig_len} -max_iteration {params.max_iteration} \
#         -thread {threads} -markerset {params.markerset}"""


#  rule get_read_counts #from jeremy zucker scripts
#     input:
#         verse_counts = rules.runcounting.output
#     output:
#         verse_read_counts = "output/{eid}/annotation/quantification/{sample}_read_counts.tsv"
#     message:
#         "Generate read frequencies from VERSE"
#     shell:
#         """python get_read_counts.py --read_count_files {input} --out {output}"""


#  rule get_function_counts #from jeremy zucker scripts
#     input:
#         read_counts = rules.getreadcounts.output
#         annotation_summary_table = rules.lcaparselast.output
#     output:
#         functional_counts_ec = "output/{eid}/annotation/quantification/{sample}_ec.tsv"
#         functional_counts_cog = "output/{eid}/annotation/quantification/{sample}_cog.tsv"
#         functional_counts_ko = "output/{eid}/annotation/quantification/{sample}_ko.tsv"
#         functional_counts_dbcan = "output/{eid}/annotation/quantification/{sample}_dbcan.tsv"
#         functional_counts_metacyc = "output/{eid}/annotation/quantification/{sample}_metacyc.tsv"
#     shell:
#         """python src/get_function_counts.py {input.annotation_summary_table} {input.read_counts} \
#            --function {wildcards.function} --out {output}"""


#  rule get_taxa_counts: #from jeremy zucker scripts
#     input:
#         annotation_summary_table = rules.lcaparselast.output
#         read_counts = rules.getreadcounts.output
#         "Taxonomy/ncbi.map",
#         "Taxonomy/nodes.dmp",
#         "Taxonomy/merged.dmp"
#     output:
#         "{sample}_Counts/{sample}_taxa_{rank}_counts.tsv"
#     shell:
#         "python src/get_rank_counts.py {input} --rank {wildcards.rank} --out {output}"


#  rule funtaxa_counts:
#     input:
#         annotation_summary_table = rules.lcaparselast.output
#         read_counts = rules.getreadcounts.output
#         "Taxonomy/ncbi.map",
#         "Taxonomy/nodes.dmp",
#         "Taxonomy/merged.dmp"
#     output:
#         "{sample}_Counts/{sample}_funtaxa_{function}_{rank}_counts.tsv"
#     shell:
#         """python src/get_fun_rank_counts.py {input} --rank {wildcards.rank} --function {wildcards.function} --out {output}"""
