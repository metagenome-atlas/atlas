import os
import re
import sys
from glob import glob
from snakemake.utils import report


def get_ribosomal_rna_input(wildcards):
    inputs = []
    data_type = config["samples"][wildcards.sample].get("type", "metagenome").lower()

    clean_reads = "{sample}/quality_control/{sample}_02_pe.fastq.gz".format(sample=wildcards.sample)
    rrna_reads = "{sample}/quality_control/{sample}_02_rRNA.fastq.gz".format(sample=wildcards.sample)

    if data_type == "metagenome" and os.path.exists(rrna_reads):
        return [clean_reads, rrna_reads]
    else:
        return [clean_reads]


def get_quality_controlled_reads(wildcards):
    # reads that have gone through ATLAS QC
    fastq = "{sample}/quality_control/{sample}_03_pe.fastq.gz".format(sample=wildcards.sample)
    # QA'd reads; the user wants to begin at assembly step
    if config.get("workflow", "complete") == "assembly":
        fastq = config["samples"][wildcards.sample]["fastq"]
    return fastq


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


def bb_cov_stats_to_maxbin(tsv_in, tsv_out):
    with open(tsv_in) as fi, open(tsv_out, "w") as fo:
        # header
        next(fi)
        for line in fi:
            toks = line.strip().split("\t")
            print(toks[0], toks[1], sep="\t", file=fo)


rule quality_filter:
    input:
        lambda wc: config["samples"][wc.sample]["fastq"]
    output:
        pe = "{sample,(?:(?!coassemblies)[\w-]+)+}/quality_control/{sample}_00_pe.fastq.gz",
        se = "{sample}/quality_control/{sample}_00_se.fastq.gz",
        stats = "{sample}/logs/{sample}_quality_filtering_stats.txt"
    benchmark:
        "logs/benchmarks/quality_filter/{sample}.txt"
    params:
        lref = "lref=%s" % config.get("preprocess_adapters") if config.get("preprocess_adapters") else "",
        rref = "rref=%s" % config.get("preprocess_adapters") if config.get("preprocess_adapters") else "",
        mink = config.get("preprocess_adapter_min_k", PREPROCESS_ADAPTER_MIN_K),
        trimq = config.get("preprocess_minimum_base_quality", PREPROCESS_MINIMUM_BASE_QUALITY),
        hdist = config.get("preprocess_allowable_kmer_mismatches", PREPROCESS_ALLOWABLE_KMER_MISMATCHES),
        k = config.get("preprocess_reference_kmer_match_length", PREPROCESS_REFERENCE_KMER_MATCH_LENGTH),
        qtrim = config.get("qtrim", QTRIM),
        minlength = config.get("preprocess_minimum_passing_read_length", PREPROCESS_MINIMUM_PASSING_READ_LENGTH),
        minbasefrequency = config.get("preprocess_minimum_base_frequency", PREPROCESS_MINIMUM_BASE_FREQUENCY),
        inputs = lambda wc: "in=%s" % config["samples"][wc.sample]["fastq"][0] if len(config["samples"][wc.sample]["fastq"]) == 1 else "in=%s in2=%s" % (config["samples"][wc.sample]["fastq"][0], config["samples"][wc.sample]["fastq"][1]),
        interleaved = lambda wc: "t" if config["samples"][wc.sample].get("paired", True) and len(config["samples"][wc.sample]["fastq"]) == 1 else "f"
    log:
        "{sample}/logs/{sample}_quality_filter.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbduk2.sh {params.inputs} out={output.pe} outs={output.se} \
               {params.rref} {params.lref} mink={params.mink} qout=33 \
               stats={output.stats} hdist={params.hdist} k={params.k} \
               trimq={params.trimq} qtrim={params.qtrim} threads={threads} \
               minlength={params.minlength} minbasefrequency={params.minbasefrequency} \
               interleaved={params.interleaved} overwrite=true 2> {log}"""


if config.get("perform_error_correction", True):
    rule error_correction:
        input:
            "{sample}/quality_control/{sample}_00_pe.fastq.gz"
        output:
            "{sample}/quality_control/{sample}_01_pe.fastq.gz"
        benchmark:
            "logs/benchmarks/error_correction/{sample}.txt"
        log:
            "{sample}/logs/{sample}_error_correction.log"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        resources:
            mem = int(re.findall(r"(\d+)", config.get("java_mem", "32"))[0])
        threads:
            config.get("threads", 1)
        shell:
            """{SHPFXM} tadpole.sh in={input} out={output} mode=correct threads={threads} \
                   ecc=t ecco=t 2> {log}"""


    # if there are no references, decontamination will be skipped
    rule decontamination:
        input:
            "{sample}/quality_control/{sample}_01_pe.fastq.gz"
        output:
            dbs = ["{sample}/quality_control/{sample}_02_%s.fastq.gz" % db for db in list(config["contaminant_references"].keys())],
            stats = "{sample}/quality_control/{sample}_decontamination_reference_stats.txt",
            clean = temp("{sample}/quality_control/{sample}_02_pe.fastq.gz")
        benchmark:
            "logs/benchmarks/decontamination/{sample}.txt"
        params:
            refs_in = " ".join(["ref_%s=%s" % (n, fa) for n, fa in config["contaminant_references"].items()]),
            refs_out = lambda wc: " ".join(["out_{ref}={sample}/quality_control/{sample}_02_{ref}.fastq.gz".format(ref=n, sample=wc.sample) for n in list(config["contaminant_references"].keys())]),
            maxindel = config.get("contaminant_max_indel", 20),
            minratio = config.get("contaminant_min_ratio", 0.65),
            minhits = config.get("contaminant_minimum_hits", 1),
            ambiguous = config.get("contaminant_ambiguous", "best"),
            k = config.get("contaminant_kmer_length", 13),
            interleaved = lambda wc: "t" if config["samples"][wc.sample].get("paired", True) else "auto"
        log:
            "{sample}/logs/{sample}_decontamination.log"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        shell:
            """{SHPFXM} bbsplit.sh nodisk=t {params.refs_in} in={input} outu={output.clean} \
                   {params.refs_out} maxindel={params.maxindel} minratio={params.minratio} \
                   minhits={params.minhits} ambiguous={params.ambiguous} refstats={output.stats}\
                   interleaved={params.interleaved} threads={threads} k={params.k} local=t 2> {log}"""


else:
    rule decontamination:
        # verify implications of no disk index when using very large reference
        input:
            "{sample}/quality_control/quality_filter/{sample}_00_pe.fastq.gz"
        output:
            dbs = ["{sample}/quality_control/{sample}_02_%s.fastq.gz" % db for db in list(config["contaminant_references"].keys())],
            stats = "{sample}/quality_control/{sample}_decontamination_reference_stats.txt",
            clean = temp("{sample}/quality_control/{sample}_02_pe.fastq.gz")
        benchmark:
            "logs/benchmarks/decontamination/{sample}.txt"
        params:
            refs_in = " ".join(["ref_%s=%s" % (n, fa) for n, fa in config["contaminant_references"].items()]),
            refs_out = lambda wc: " ".join(["out_{ref}={sample}/quality_control/{sample}_{ref}.fastq.gz".format(ref=n, sample=wc.sample) for n in list(config["contaminant_references"].keys())]),
            maxindel = config.get("contaminant_max_indel", 20),
            minratio = config.get("contaminant_min_ratio", 0.65),
            minhits = config.get("contaminant_minimum_hits", 1),
            ambiguous = config.get("contaminant_ambiguous", "best"),
            k = config.get("contaminant_kmer_length", 13),
            interleaved = lambda wc: "t" if config["samples"][wc.sample].get("paired", True) else "auto"
        log:
            "{sample}/logs/{sample}_decontamination.log"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        shell:
            """{SHPFXM} bbsplit.sh nodisk=t {params.refs_in} in={input} outu={output.clean} \
                   {params.refs_out} maxindel={params.maxindel} minratio={params.minratio} \
                   minhits={params.minhits} ambiguous={params.ambiguous} refstats={output.stats}\
                   interleaved={params.interleaved} threads={threads} k={params.k} local=t 2> {log}"""


rule postprocess_after_decontamination:
    input:
        get_ribosomal_rna_input
    output:
        "{sample}/quality_control/{sample}_03_pe.fastq.gz"
    threads:
        1
    shell:
        "{SHPFXS} cat {input} > {output}"


rule normalize_coverage_across_kmers:
    input:
        get_quality_controlled_reads
    output:
        "{sample}/quality_control/{sample}_04_%s_pe.fastq.gz" % NORMALIZATION
    benchmark:
        "logs/benchmarks/normalization/{sample}.txt"
    params:
        k = config.get("normalization_kmer_length", NORMALIZATION_KMER_LENGTH),
        t = config.get("normalization_target_depth", NORMALIZATION_TARGET_DEPTH),
        minkmers = config.get("normalization_minimum_kmers", NORMALIZATION_MINIMUM_KMERS),
        inputs = get_bbtools_input, lambda wc: "in=%s" % config["samples"][wc.sample]["fastq"][0] if len(config["samples"][wc.sample]["fastq"]) == 1 else "in=%s in2=%s" % (config["samples"][wc.sample]["fastq"][0], config["samples"][wc.sample]["fastq"][1]),
        interleaved = lambda wc: "t" if config["samples"][wc.sample].get("paired", True) and len(config["samples"][wc.sample]["fastq"]) == 1 else "f"
    log:
        "{sample}/logs/{sample}_%s.log" % NORMALIZATION
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbnorm.sh in={input} out={output} k={params.k} t={params.t} \
               interleaved={params.interleaved} minkmers={params.minkmers} prefilter=t \
               threads={threads} 2> {log}"""


if config.get("assembler", "megahit") == "megahit":
    rule run_megahit:
        input:
            "{sample}/quality_control/{sample}_04_%s_pe.fastq.gz" % NORMALIZATION
        output:
            temp("{sample}/{assembler}/{sample}_prefilter.contigs.fa")
        benchmark:
            "logs/benchmarks/{assembler}/{sample}.txt"
        shadow:
            "full"
        params:
            read_flag = lambda wc: "--12" if config["samples"][wc.sample].get("paired", True) else "--read",
            memory = config.get("megahit_memory", 0.90),
            min_count = config.get("megahit_min_count", 2),
            k_min = config.get("megahit_k_min", 21),
            k_max = config.get("megahit_k_max", 121),
            k_step = config.get("megahit_k_step", 20),
            merge_level = config.get("megahit_merge_level", "20,0.98"),
            prune_level = config.get("megahit_prune_level", 2),
            low_local_ratio = config.get("megahit_low_local_ratio", 0.2),
            min_contig_len = config.get("minimum_contig_length", 1000),
            outdir = lambda wc: "{sample}/{assembler}".format(sample=wc.sample, assembler=ASSEMBLER)
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        shell:
            """{SHPFXM} megahit --continue --tmp-dir {TMPDIR} --num-cpu-threads {threads} {params.read_flag} {input} \
                   --k-min {params.k_min} --k-max {params.k_max} --k-step {params.k_step} \
                   --out-dir {params.outdir} --out-prefix {wildcards.sample}_prefilter \
                   --min-contig-len {params.min_contig_len} --min-count {params.min_count} \
                   --merge-level {params.merge_level} --prune-level {params.prune_level} \
                   --low-local-ratio {params.low_local_ratio}"""


    rule rename_megahit_output:
        input:
            "{sample}/{assembler}/{sample}_prefilter.contigs.fa"
        output:
            temp("{sample}/{assembler}/{sample}_raw_contigs.fasta")
        shell:
            "{SHPFXS} cp {input} {output}"

else:
    rule run_spades:
        input:
            "{sample}/quality_control/{sample}_04_%s_pe.fastq.gz" % NORMALIZATION
        output:
            temp("{sample}/{assembler}/contigs.fasta")
        benchmark:
            "logs/benchmarks/{assembler}/{sample}.txt"
        params:
            # memory = config["assembly"].get("memory", 0.90)
            read_flag = lambda wc: "--12" if config["samples"][wc.sample].get("paired", True) else "-s",
            k = config.get("spades_k", "auto"),
            outdir = lambda wc: "{sample}/{assembler}".format(sample=wc.sample, assembler=ASSEMBLER)
        log:
            "{sample}/{assembler}/spades.log"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        shell:
            """{SHPFXM} spades.py -t {threads} -o {params.outdir} --meta {params.read_flag} {input}"""


    rule rename_spades_output:
        input:
            "{sample}/{assembler}/contigs.fasta"
        output:
            temp("{sample}/{assembler}/{sample}_raw_contigs.fasta")
        shell:
            "{SHPFXS} cp {input} {output}"


rule rename_contigs:
    # standardizes header labels within contig FASTAs
    input:
        "{sample}/{assembler}/{sample}_raw_contigs.fasta"
    output:
        "{sample}/{assembler}/{sample}_prefilter_contigs.fasta"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    shell:
        """rename.sh in={input} out={output} ow=t prefix={wildcards.sample}"""


rule calculate_prefiltered_contigs_stats:
    input:
        "{sample}/{assembler}/{sample}_prefilter_contigs.fasta"
    output:
        "{sample}/{assembler}/contig_stats/prefilter_contig_stats.txt"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        1
    shell:
        "{SHPFXS} stats.sh in={input} format=3 > {output}"


rule calculate_prefiltered_contig_coverage_stats:
    input:
        fasta = "{sample}/{assembler}/{sample}_prefilter_contigs.fasta",
        fastq = get_quality_controlled_reads
    output:
        bhist = "{sample}/{assembler}/contig_stats/prefilter_base_composition.txt",
        bqhist = "{sample}/{assembler}/contig_stats/prefilter_box_quality.txt",
        mhist = "{sample}/{assembler}/contig_stats/prefilter_mutation_rates.txt",
        statsfile = "{sample}/{assembler}/contig_stats/prefilter_mapping_stats.txt",
        covstats = "{sample}/{assembler}/contig_stats/prefilter_coverage_stats.txt"
    benchmark:
        "logs/benchmarks/calculate_prefiltered_contig_coverage_stats/{sample}.txt"
    params:
        interleaved = lambda wc: "t" if config["samples"][wc.sample].get("paired", True) else "auto"
    log:
        "{sample}/{assembler}/logs/prefiltered_contig_coverage_stats.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbmap.sh nodisk=t ref={input.fasta} in={input.fastq} fast=t \
               interleaved={params.interleaved} threads={threads} bhist={output.bhist} \
               bqhist={output.bqhist} mhist={output.mhist} statsfile={output.statsfile} \
               covstats={output.covstats} 2> {log}"""


rule filter_by_coverage:
    input:
        fasta = "{sample}/{assembler}/{sample}_prefilter_contigs.fasta",
        covstats = "{sample}/{assembler}/contig_stats/prefilter_coverage_stats.txt"
    output:
        fasta = "{sample}/{assembler}/{sample}_contigs.fasta",
        removed_names = "{sample}/{assembler}/{sample}_discarded_contigs.fasta"
    params:
        minc = config.get("minimum_average_coverage", 5),
        minp = config.get("minimum_percent_covered_bases", 40),
        minr = config.get("minimum_mapped_reads", 0),
        minl = config.get("minimum_contig_length", 1000),
        trim = config.get("contig_trim_bp", 0)
    log:
        "{sample}/{assembler}/logs/filter_by_coverage.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        1
    shell:
        """{SHPFXS} filterbycoverage.sh in={input.fasta} cov={input.covstats} out={output.fasta} \
               outd={output.removed_names} minc={params.minc} minp={params.minp} \
               minr={params.minr} minl={params.minl} trim={params.trim} 2> {log}"""


rule align_reads_to_filtered_contigs:
    input:
        fasta = "{sample}/{assembler}/{sample}_contigs.fasta",
        fastq = get_quality_controlled_reads
    output:
        sam = temp("{sample}/{assembler}/alignments/{sample}.sam"),
        bhist = "{sample}/{assembler}/contig_stats/postfilter_base_composition.txt",
        bqhist = "{sample}/{assembler}/contig_stats/postfilter_box_quality.txt",
        mhist = "{sample}/{assembler}/contig_stats/postfilter_mutation_rates.txt",
        gchist = "{sample}/{assembler}/contig_stats/postfilter_gc_rates.txt",
        statsfile = "{sample}/{assembler}/contig_stats/postfilter_mapping_stats.txt",
        covstats = "{sample}/{assembler}/contig_stats/postfilter_coverage_stats.txt"
    benchmark:
        "logs/benchmarks/align_reads_to_filtered_contigs/{sample}.txt"
    params:
        interleaved = lambda wc: "t" if config["samples"][wc.sample].get("paired", True) else "auto",
        maxsites = config.get("maximum_counted_map_sites", 10)
    log:
        "{sample}/{assembler}/logs/contig_coverage_stats.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbmap.sh nodisk=t ref={input.fasta} in={input.fastq} trimreaddescriptions=t \
               out={output.sam} mappedonly=t threads={threads} bhist={output.bhist} \
               bqhist={output.bqhist} mhist={output.mhist} gchist={output.gchist} \
               statsfile={output.statsfile} covstats={output.covstats} mdtag=t xstag=fs nmtag=t \
               sam=1.3 local=t ambiguous=all interleaved={params.interleaved} secondary=t ssao=t \
               maxsites={params.maxsites} 2> {log}"""


if config.get("perform_genome_binning", True):
    rule make_maxbin_abundance_file:
        input:
            "{sample}/{assembler}/contig_stats/postfilter_coverage_stats.txt"
        output:
            "{sample}/{assembler}/genomic_bins/{sample}_contig_coverage.tsv"
        run:
            bb_cov_stats_to_maxbin(input[0], output[0])


    rule run_maxbin:
        input:
            fasta = "{sample}/{assembler}/{sample}_contigs.fasta",
            abundance = "{sample}/{assembler}/genomic_bins/{sample}_contig_coverage.tsv"
        output:
            # fastas will need to be dynamic if we do something with them at a later time
            summary = "{sample}/{assembler}/genomic_bins/{sample}.summary",
            marker = "{sample}/{assembler}/genomic_bins/{sample}.marker"
        benchmark:
            "logs/benchmarks/maxbin2/{sample}.txt"
        params:
            mi = config.get("maxbin_max_iteration", 50),
            mcl = config.get("maxbin_min_contig_length", 200),
            pt = config.get("maxbin_prob_threshold", 0.9)
        log:
            "{sample}/{assembler}/logs/maxbin2.log"
        conda:
            "%s/optional_genome_binning.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        shell:
            """{SHPFXM} run_MaxBin.pl -contig {input.fasta} -abund {input.abundance} \
                -out {wildcards.sample}/{wildcards.assembler}/genomic_bins/{wildcards.sample} \
                -min_contig_length {params.mcl} -thread {threads} -prob_threshold {params.pt} \
                -max_iteration {params.mi} > {log}"""


    rule initialize_checkm:
        # input:
        output:
            "%s/test_data/637000110.fna" % CHECKMDIR,
            "%s/taxon_marker_sets.tsv" % CHECKMDIR,
            "%s/selected_marker_sets.tsv" % CHECKMDIR,
            "%s/pfam/tigrfam2pfam.tsv" % CHECKMDIR,
            "%s/pfam/Pfam-A.hmm.dat" % CHECKMDIR,
            "%s/img/img_metadata.tsv" % CHECKMDIR,
            "%s/hmms_ssu/SSU_euk.hmm" % CHECKMDIR,
            "%s/hmms_ssu/SSU_bacteria.hmm" % CHECKMDIR,
            "%s/hmms_ssu/SSU_archaea.hmm" % CHECKMDIR,
            "%s/hmms_ssu/createHMMs.py" % CHECKMDIR,
            "%s/hmms/phylo.hmm.ssi" % CHECKMDIR,
            "%s/hmms/phylo.hmm" % CHECKMDIR,
            "%s/hmms/checkm.hmm.ssi" % CHECKMDIR,
            "%s/hmms/checkm.hmm" % CHECKMDIR,
            "%s/genome_tree/missing_duplicate_genes_97.tsv" % CHECKMDIR,
            "%s/genome_tree/missing_duplicate_genes_50.tsv" % CHECKMDIR,
            "%s/genome_tree/genome_tree.taxonomy.tsv" % CHECKMDIR,
            "%s/genome_tree/genome_tree_reduced.refpkg/phylo_modelJqWx6_.json" % CHECKMDIR,
            "%s/genome_tree/genome_tree_reduced.refpkg/genome_tree.tre" % CHECKMDIR,
            "%s/genome_tree/genome_tree_reduced.refpkg/genome_tree.log" % CHECKMDIR,
            "%s/genome_tree/genome_tree_reduced.refpkg/genome_tree.fasta" % CHECKMDIR,
            "%s/genome_tree/genome_tree_reduced.refpkg/CONTENTS.json" % CHECKMDIR,
            "%s/genome_tree/genome_tree.metadata.tsv" % CHECKMDIR,
            "%s/genome_tree/genome_tree_full.refpkg/phylo_modelEcOyPk.json" % CHECKMDIR,
            "%s/genome_tree/genome_tree_full.refpkg/genome_tree.tre" % CHECKMDIR,
            "%s/genome_tree/genome_tree_full.refpkg/genome_tree.log" % CHECKMDIR,
            "%s/genome_tree/genome_tree_full.refpkg/genome_tree.fasta" % CHECKMDIR,
            "%s/genome_tree/genome_tree_full.refpkg/CONTENTS.json" % CHECKMDIR,
            "%s/genome_tree/genome_tree.derep.txt" % CHECKMDIR,
            "%s/.dmanifest" % CHECKMDIR,
            "%s/distributions/td_dist.txt" % CHECKMDIR,
            "%s/distributions/gc_dist.txt" % CHECKMDIR,
            "%s/distributions/cd_dist.txt" % CHECKMDIR,
            touched_output = "logs/checkm_init.txt"
        params:
            database_dir = CHECKMDIR
        conda:
            "%s/optional_genome_binning.yaml" % CONDAENV
        script:
            "initialize_checkm.py"


    rule run_checkm_lineage_wf:
        input:
            init_checkm = "%s/hmms/checkm.hmm" % CHECKMDIR,
            bins = "{sample}/{assembler}/genomic_bins/{sample}.marker"
        output:
            "{sample}/{assembler}/genomic_bins/checkm/completeness.tsv"
        params:
            bin_dir = "{sample}/{assembler}/genomic_bins",
            output_dir = "{sample}/{assembler}/genomic_bins/checkm"
        conda:
            "%s/optional_genome_binning.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        shell:
            """rm -r {params.output_dir} && {SHPFXM} checkm lineage_wf --file {params.output_dir}/completeness.tsv --tab_table \
                   --quiet --extension fasta --threads {threads} {params.bin_dir} \
                   {params.output_dir}"""


    rule run_checkm_tree_qa:
        input:
            "{sample}/{assembler}/genomic_bins/checkm/completeness.tsv"
        output:
            "{sample}/{assembler}/genomic_bins/checkm/taxonomy.tsv"
        params:
            output_dir = "{sample}/{assembler}/genomic_bins/checkm"
        conda:
            "%s/optional_genome_binning.yaml" % CONDAENV
        shell:
            """{SHPFXS} checkm tree_qa --tab_table --out_format 2 \
                   --file {params.output_dir}/taxonomy.tsv {params.output_dir}"""


    # rule compile_bin_tax_assignments:
    #     input:
    #         "{sample}/{assembler}/genomic_bins/checkm/completeness.tsv",
    #         "{sample}/{assembler}/genomic_bins/checkm/taxonomy.tsv"
    #     output:
    #         "{sample}/{assembler}/genomic_bins/checkm/completeness_and_taxonomy.tsv"



rule convert_sam_to_bam:
    input:
        "{sample}/{assembler}/alignments/{sample}.sam"
    output:
        temp("{sample}/{assembler}/alignments/{sample}.bam")
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} samtools view -@ {threads} -bSh1 {input} \
               | samtools sort -m 1536M -@ {threads} -T {TMPDIR}/{wildcards.sample}_tmp -o {output} -O bam -"""


rule create_bam_index:
    input:
        "{sample}/{assembler}/alignments/{sample}.bam"
    output:
        temp("{sample}/{assembler}/alignments/{sample}.bam.bai")
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        1
    shell:
        "{SHPFXS} samtools index {input}"


rule calculate_final_contigs_stats:
    input:
        "{sample}/{assembler}/{sample}_contigs.fasta"
    output:
        "{sample}/{assembler}/contig_stats/final_contig_stats.txt"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        1
    shell:
        "{SHPFXS} stats.sh in={input} format=3 > {output}"


rule run_prokka_annotation:
    input:
        "{sample}/{assembler}/{sample}_contigs.fasta"
    output:
        discrepancy = "{sample}/{assembler}/functional_annotation/prokka/{sample}.err",
        faa = "{sample}/{assembler}/functional_annotation/prokka/{sample}.faa",
        ffn = "{sample}/{assembler}/functional_annotation/prokka/{sample}.ffn",
        fna = "{sample}/{assembler}/functional_annotation/prokka/{sample}.fna",
        fsa = "{sample}/{assembler}/functional_annotation/prokka/{sample}.fsa",
        gbk = "{sample}/{assembler}/functional_annotation/prokka/{sample}.gbk",
        gff = "{sample}/{assembler}/functional_annotation/prokka/{sample}.gff",
        log = "{sample}/{assembler}/functional_annotation/prokka/{sample}.log",
        sqn = "{sample}/{assembler}/functional_annotation/prokka/{sample}.sqn",
        tbl = "{sample}/{assembler}/functional_annotation/prokka/{sample}.tbl",
        tsv = "{sample}/{assembler}/functional_annotation/prokka/{sample}.tsv",
        txt = "{sample}/{assembler}/functional_annotation/prokka/{sample}.txt"
    benchmark:
        "logs/benchmarks/prokka/{sample}.txt"
    params:
        outdir = "{sample}/{assembler}/functional_annotation/prokka",
        kingdom = config.get("prokka_kingdom", "Bacteria")
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} prokka --outdir {params.outdir} --force --prefix {wildcards.sample} \
               --locustag {wildcards.sample} --kingdom {params.kingdom} --metagenome \
               --cpus {threads} {input}"""


rule update_prokka_tsv:
    input:
        "{sample}/{assembler}/functional_annotation/prokka/{sample}.gff"
    output:
        "{sample}/{assembler}/functional_annotation/prokka/{sample}_plus.tsv"
    shell:
        """{SHPFXS} atlas gff2tsv {input} {output}"""


rule convert_gff_to_gtf:
    input:
        "{sample}/{assembler}/functional_annotation/prokka/{sample}.gff"
    output:
        "{sample}/{assembler}/functional_annotation/prokka/{sample}.gtf"
    run:
        gff_to_gtf(input[0], output[0])


rule remove_pcr_duplicates:
    input:
        bam = "{sample}/{assembler}/alignments/{sample}.bam",
        bai = "{sample}/{assembler}/alignments/{sample}.bam.bai"
    output:
        bam = "{sample}/{assembler}/alignments/{sample}_markdup.bam",
        txt = "{sample}/{assembler}/alignments/{sample}_markdup_metrics.txt"
    benchmark:
        "logs/benchmarks/picard_mark_duplicates/{sample}.txt"
    params:
        java_mem = config.get("java_mem", "32g")
    conda:
        "%s/required_packages.yaml" % CONDAENV
    resources:
        mem = int(re.findall(r"(\d+)", config.get("java_mem", "32"))[0])
    shell:
        """{SHPFXS} picard MarkDuplicates -Xmx{params.java_mem} INPUT={input.bam} \
               OUTPUT={output.bam} METRICS_FILE={output.txt} ASSUME_SORT_ORDER=coordinate \
               MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=TRUE \
               VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE"""


rule find_counts_per_region:
    input:
        gtf = "{sample}/{assembler}/functional_annotation/prokka/{sample}.gtf",
        bam = "{sample}/{assembler}/alignments/{sample}_markdup.bam"
    output:
        summary = "{sample}/{assembler}/functional_annotation/feature_counts/{sample}_counts.txt.summary",
        counts = "{sample}/{assembler}/functional_annotation/feature_counts/{sample}_counts.txt"
    params:
        min_read_overlap = config.get("minimum_region_overlap", 1),
        paired_mode = lambda wc: "-p" if config["samples"][wc.sample].get("paired", True) else "",
        multi_mapping = "-M" if config.get("count_multi_mapped_reads") else "",
        primary_only = "--primary" if config.get("primary_only", False) else ""
    log:
        "{sample}/{assembler}/logs/counts_per_region.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} featureCounts {params.paired_mode} -F gtf -T {threads} \
               {params.multi_mapping} -t CDS -g ID -a {input.gtf} -o {output.counts} \
               {input.bam} 2> {log}"""


rule run_diamond_blastp:
    input:
        fasta = "{sample}/{assembler}/functional_annotation/prokka/{sample}.faa",
        db = config["diamond_db"]
    output:
        "{sample}/{assembler}/functional_annotation/refseq/{sample}_hits.tsv"
    benchmark:
        "logs/benchmarks/run_diamond_blastp/{sample}.txt"
    params:
        tmpdir = "--tmpdir %s" % TMPDIR if TMPDIR else "",
        top_seqs = config.get("diamond_top_seqs", 2),
        e_value = config.get("diamond_e_value", 0.000001),
        min_identity = config.get("diamond_min_identity", 50),
        query_cover = config.get("diamond_query_coverage", 60),
        gap_open = config.get("diamond_gap_open", 11),
        gap_extend = config.get("diamond_gap_extend", 1),
        block_size = config.get("diamond_block_size", 2),
        index_chunks = config.get("diamond_index_chunks", 4),
        run_mode = "--more-sensitive" if not config.get("diamond_run_mode", "") == "fast" else ""
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} diamond blastp --threads {threads} --outfmt 6 --out {output} \
               --query {input.fasta} --db {input.db} --top {params.top_seqs} \
               --evalue {params.e_value} --id {params.min_identity} \
               --query-cover {params.query_cover} {params.run_mode} --gapopen {params.gap_open} \
               --gapextend {params.gap_extend} {params.tmpdir} --block-size {params.block_size} \
               --index-chunks {params.index_chunks}"""


rule add_contig_metadata:
    input:
        hits = "{sample}/{assembler}/functional_annotation/refseq/{sample}_hits.tsv",
        gff = "{sample}/{assembler}/functional_annotation/prokka/{sample}.gff"
    output:
        temp("{sample}/{assembler}/functional_annotation/refseq/{sample}_hits_plus.tsv")
    shell:
        "{SHPFXS} atlas munge-blast {input.hits} {input.gff} {output}"


rule sort_munged_blast_hits:
    # ensure blast hits are grouped by contig, ORF, and then decreasing by bitscore
    input:
        "{sample}/{assembler}/functional_annotation/refseq/{sample}_hits_plus.tsv"
    output:
        "{sample}/{assembler}/functional_annotation/refseq/{sample}_hits_plus_sorted.tsv"
    shell:
        "{SHPFXS} sort -k1,1 -k2,2 -k13,13rn {input} > {output}"


rule parse_blastp:
    # assign a taxonomy to contigs using the consensus of the ORF assignments
    input:
        "{sample}/{assembler}/functional_annotation/refseq/{sample}_hits_plus_sorted.tsv"
    output:
        "{sample}/{assembler}/functional_annotation/refseq/{sample}_tax_assignments.tsv"
    params:
        namemap = config["refseq_namemap"],
        treefile = config["refseq_tree"],
        summary_method = config.get("summary_method", "lca"),
        aggregation_method = config.get("aggregation_method", "lca-majority"),
        majority_threshold = config.get("majority_threshold", 0.51),
        min_identity = config.get("diamond_min_identity", 50),
        min_bitscore = config.get("min_bitscore", 0),
        min_length = config.get("min_length", 20),
        max_evalue = config.get("diamond_e_value", 0.000001),
        max_hits = config.get("max_hits", 100),
        top_fraction = (100 - config.get("diamond_top_seqs", 5)) * 0.01
    shell:
        """{SHPFXS} atlas refseq --summary-method {params.summary_method} \
               --aggregation-method {params.aggregation_method} \
               --majority-threshold {params.majority_threshold} \
               --min-identity {params.min_identity} \
               --min-bitscore {params.min_bitscore} \
               --min-length {params.min_length} \
               --max-evalue {params.max_evalue} \
               --max-hits {params.max_hits} \
               --top-fraction {params.top_fraction} \
               {input} {params.namemap} {params.treefile} {output}"""


if config.get("perform_genome_binning", True):
    rule merge_sample_tables:
        input:
            prokka = "{sample}/{assembler}/functional_annotation/prokka/{sample}_plus.tsv",
            refseq = "{sample}/{assembler}/functional_annotation/refseq/{sample}_tax_assignments.tsv",
            counts = "{sample}/{assembler}/functional_annotation/feature_counts/{sample}_counts.txt",
            completeness = "{sample}/{assembler}/genomic_bins/checkm/completeness.tsv",
            taxonomy = "{sample}/{assembler}/genomic_bins/checkm/taxonomy.tsv"
        output:
            "{sample}/{assembler}/{sample}_annotations.txt"
        params:
            fastas = lambda wc: " --fasta ".join(glob("{sample}/{assembler}/genomic_bins/{sample}.*.fasta".format(sample=wc.sample, assembler=wc.assembler)))
        shell:
            "{SHPFXS} atlas merge-tables --counts {input.counts} \
                 --completeness {input.completeness} --taxonomy {input.taxonomy} \
                 --fasta {params.fastas} {input.prokka} {input.refseq} {output}"


else:
    rule merge_sample_tables:
        input:
            prokka = "{sample}/{assembler}/functional_annotation/prokka/{sample}_plus.tsv",
            refseq = "{sample}/{assembler}/functional_annotation/refseq/{sample}_tax_assignments.tsv",
            counts = "{sample}/{assembler}/functional_annotation/feature_counts/{sample}_counts.txt"
        output:
            "{sample}/{assembler}/{sample}_annotations.txt"
        shell:
            "{SHPFXS} atlas merge-tables --counts {input.counts} {input.prokka} {input.refseq} {output}"


# rule assembly_report:
#
#     input:
#         contig_stats = "{sample}/{assembler}/contig_stats/final_contig_stats.txt",
#         base_comp = "{sample}/{assembler}/contig_stats/postfilter_base_composition.txt"
#         # css = os.path.join(workflow.basedir, "resources", "report.css")
#     output:
#         html = "{sample}/{sample}_{assembler}_README.html"
#     shadow:
#         "shallow"
#     run:
#         import pandas as pd
#         # contig stats table
#         df = pd.read_csv(input.contig_stats, sep="\t")
#         contig_stats_csv = "contig_stats.csv"
#         df.to_csv(contig_stats_csv,
#                   columns=["n_contigs", "contig_bp", "ctg_N50", "ctg_N90", "ctg_max", "gc_avg"],
#                   index=False)
#
#         # read base composition across final contigs
#         df = pd.read_csv(input.base_comp, sep="\t")
#         base_composition_positions = "['%s']" % "', '".join(map(str, df["#Pos"]))
#         base_composition_a = "[%s]" % ", ".join(map(str, df["A"]))
#         base_composition_c = "[%s]" % ", ".join(map(str, df["C"]))
#         base_composition_g = "[%s]" % ", ".join(map(str, df["G"]))
#         base_composition_t = "[%s]" % ", ".join(map(str, df["T"]))
#         base_composition_n = "[%s]" % ", ".join(map(str, df["N"]))
#
#         report("""
#
# ===========================================================================================
# Sample Report - Sample: {wildcards.sample}
# ===========================================================================================
#
# .. raw:: html
#
#     body{font-family:Helvetica,arial,sans-serif;font-size:14px;line-height:1.6;background-color:#fff;padding:30px;color:#333}body > :first-child{margin-top:0!important}body > :last-child{margin-bottom:0!important}a{color:#4183C4;text-decoration:none}a.absent{color:#c00}a.anchor{display:block;padding-left:30px;margin-left:-30px;cursor:pointer;position:absolute;top:0;left:0;bottom:0}h1,h2,h3,h4,h5,h6{margin:20px 0 10px;padding:0;font-weight:700;-webkit-font-smoothing:antialiased;cursor:text;position:relative}h2:first-child,h1:first-child,h1:first-child + h2,h3:first-child,h4:first-child,h5:first-child,h6:first-child{margin-top:0;padding-top:0}h1:hover a.anchor,h2:hover a.anchor,h3:hover a.anchor,h4:hover a.anchor,h5:hover a.anchor,h6:hover a.anchor{text-decoration:none}h1 tt,h1 code{font-size:inherit}h2 tt,h2 code{font-size:inherit}h3 tt,h3 code{font-size:inherit}h4 tt,h4 code{font-size:inherit}h5 tt,h5 code{font-size:inherit}h6 tt,h6 code{font-size:inherit}h1{font-size:28px;color:#000}h2{font-size:24px;border-bottom:1px solid #ccc;color:#000}h3{font-size:18px}h4{font-size:16px}h5{font-size:14px}h6{color:#777;font-size:14px}p,blockquote,ul,ol,dl,li,table,pre{margin:15px 0}hr{background:transparent url(http://tinyurl.com/bq5kskr) repeat-x 0 0;border:0 none;color:#ccc;height:4px;padding:0}body > h2:first-child{margin-top:0;padding-top:0}body > h1:first-child{margin-top:0;padding-top:0}body > h1:first-child + h2{margin-top:0;padding-top:0}body > h3:first-child,body > h4:first-child,body > h5:first-child,body > h6:first-child{margin-top:0;padding-top:0}a:first-child h1,a:first-child h2,a:first-child h3,a:first-child h4,a:first-child h5,a:first-child h6{margin-top:0;padding-top:0}h1 p,h2 p,h3 p,h4 p,h5 p,h6 p{margin-top:0}li p.first{display:inline-block}ul,ol{padding-left:30px}ul :first-child,ol :first-child{margin-top:0}ul :last-child,ol :last-child{margin-bottom:0}dl{padding:0}dl dt{font-size:14px;font-weight:700;font-style:italic;padding:0;margin:15px 0 5px}dl dt:first-child{padding:0}dl dt > :first-child{margin-top:0}dl dt > :last-child{margin-bottom:0}dl dd{margin:0 0 15px;padding:0 15px}dl dd > :first-child{margin-top:0}dl dd > :last-child{margin-bottom:0}blockquote{border-left:4px solid #ddd;padding:0 15px;color:#777}blockquote > :first-child{margin-top:0}blockquote > :last-child{margin-bottom:0}table{padding:0;border-spacing:0;border-collapse:collapse}table tr{border-top:1px solid #ccc;background-color:#fff;margin:0;padding:0}table tr:nth-child(2n){background-color:#f8f8f8}table tr th{font-weight:700;border:1px solid #ccc;text-align:left;margin:0;padding:6px 13px}table tr td{border:1px solid #ccc;text-align:left;margin:0;padding:6px 13px}table tr th :first-child,table tr td :first-child{margin-top:0}table tr th :last-child,table tr td :last-child{margin-bottom:0}img{max-width:100%}span.frame{display:block;overflow:hidden}span.frame > span{border:1px solid #ddd;display:block;float:left;overflow:hidden;margin:13px 0 0;padding:7px;width:auto}span.frame span img{display:block;float:left}span.frame span span{clear:both;color:#333;display:block;padding:5px 0 0}span.align-center{display:block;overflow:hidden;clear:both}span.align-center > span{display:block;overflow:hidden;margin:13px auto 0;text-align:center}span.align-center span img{margin:0 auto;text-align:center}span.align-right{display:block;overflow:hidden;clear:both}span.align-right > span{display:block;overflow:hidden;margin:13px 0 0;text-align:right}span.align-right span img{margin:0;text-align:right}span.float-left{display:block;margin-right:13px;overflow:hidden;float:left}span.float-left span{margin:13px 0 0}span.float-right{display:block;margin-left:13px;overflow:hidden;float:right}span.float-right > span{display:block;overflow:hidden;margin:13px auto 0;text-align:right}code,tt{margin:0 2px;padding:0 5px;white-space:nowrap;border:1px solid #eaeaea;background-color:#f8f8f8;border-radius:3px}pre code{margin:0;padding:0;white-space:pre;border:none;background:transparent}.highlight pre{background-color:#f8f8f8;border:1px solid #ccc;font-size:13px;line-height:19px;overflow:auto;padding:6px 10px;border-radius:3px}pre{background-color:#f8f8f8;border:1px solid #ccc;font-size:13px;line-height:19px;overflow:auto;padding:6px 10px;border-radius:3px}pre code,pre tt{background-color:transparent;border:none}div#metadata{text-align:right}
#
#     <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>
#     <script src="https://code.highcharts.com/highcharts.js"></script>
#     <script src="https://code.highcharts.com/modules/exporting.js"></script>
#     <script type="text/javascript">
#
#     $(function () {{
#         $('#read_composition').highcharts({{
#             title: {{text: 'Read Base Composition by Position'}},
#             xAxis: {{title: {{text: "Position"}}, categories: {base_composition_positions}}},
#             yAxis: {{min: 0, title: {{text: 'Fraction'}}}},
#             tooltip: {{}},
#             credits: {{enabled: false}},
#             legend: {{layout: 'vertical', align: 'right', verticalAlign: 'middle', borderWidth: 0}},
#             plotOptions: {{series: {{ marker: {{ enabled: false }} }}, column: {{pointPadding: 0.2, borderWidth: 0}}}},
#             series: [{{name: 'A', data: {base_composition_a}}},
#                      {{name: 'C', data: {base_composition_c}}},
#                      {{name: 'G', data: {base_composition_g}}},
#                      {{name: 'T', data: {base_composition_t}}},
#                      {{name: 'N', data: {base_composition_n}}}]
#             }});
#     }});
#     </script>
#
# .. contents:: Contents
#     :backlinks: none
#
# Read Summary
# ------------
#
# .. raw:: html
#
#     <div id="read_composition" style="min-width: 310px; height: 500px; margin: 0 auto"></div>
#
#
# Contig Summary
# --------------
#
# .. csv-table::
#     :header-rows: 1
#     :file: {contig_stats_csv}
#
#
#                """, output.html, metadata="Author: " + config.get("author", "ATLAS"),
#                stylesheet=None, contig_stats=input.contig_stats)
