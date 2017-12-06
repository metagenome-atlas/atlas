import os
import re
import sys
from glob import glob
from snakemake.utils import report
import warnings


localrules: rename_megahit_output,rename_spades_output,initialize_checkm


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

            

# define which steps are defined before assembly and in which order

assembly_preprocessing_steps=['normalized','errorcorr','merged']

if not (config.get('merge_pairs_before_assembly',True) and paired_end ):
    assembly_preprocessing_steps.remove('merged')
assembly_preprocessing_steps=".".join(assembly_preprocessing_steps)





rule normalize_coverage_across_kmers:
    input:
        unpack(get_quality_controlled_reads) #expect SE or R1,R2 or R1,R2,SE
    output:
        temp(expand("{{sample}}/assembly/reads/normalized_{fraction}.fastq.gz",
                fraction= multifile_fractions))
    benchmark:
        "logs/benchmarks/normalization/{sample}.txt"
    params:
        k = config.get("normalization_kmer_length", NORMALIZATION_KMER_LENGTH),
        t = config.get("normalization_target_depth", NORMALIZATION_TARGET_DEPTH),
        minkmers = config.get("normalization_minimum_kmers", NORMALIZATION_MINIMUM_KMERS),
        input_single = lambda wc, input: "in=%s" % input.se if hasattr(input,'se') else "null",
        extra_single = lambda wc, input: "extra=%s,%s" % (input.R1, input.R2) if hasattr(input,'R1') else "",
        has_paired_end_files= lambda wc, input: "t" if hasattr(input,'R1') else "f",
        input_paired = lambda wc, input: "in=%s in2=%s" % (input.R1, input.R2) if hasattr(input,'R1') else "null",
        extra_paired = lambda wc, input: "extra=%s" % input.se if hasattr(input,'se') else "",
        output_single = lambda wc,output,input: "out=%s" % output[2] if hasattr(input,'R1') else "out=%s" % output[0],
        output_paired = lambda wc,output,input: "out=%s out2=%s" % (output[0],output[1]) if hasattr(input,'R1') else "null",
        interleaved = "f" #lambda wc, input: "t" if (wc.fraction=='pe') else "f"   # I don't know how to handle interleaved files at this stage
    log:
        "{sample}/logs/{sample}_normalization.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", JAVA_MEM)
    shell:
        """
            if [ {params.input_single} != "null" ];
            then
        {SHPFXM} bbnorm.sh {params.input_single} \
                {params.extra_single} \
                {params.output_single} \
                k={params.k} t={params.t} \
                interleaved={params.interleaved} minkmers={params.minkmers} prefilter=t \
                threads={threads} \
                -Xmx{resources.mem}G 2> {log}
            fi


            if [ {params.has_paired_end_files} = "t" ];
            then
        {SHPFXM} bbnorm.sh {params.input_paired} \
                {params.extra_paired} \
                {params.output_paired} \
                k={params.k} t={params.t} \
                interleaved={params.interleaved} minkmers={params.minkmers} prefilter=t \
                threads={threads} \
                -Xmx{resources.mem}G 2>> {log}
            fi

            """

rule error_correction:
    input:
        expand("{{sample}}/assembly/reads/{{previous_steps}}_{fraction}.fastq.gz",
            fraction=multifile_fractions)
    output:
        temp(expand("{{sample}}/assembly/reads/{{previous_steps}}.errorcorr_{fraction}.fastq.gz",
            fraction=multifile_fractions))
    benchmark:
        "logs/benchmarks/error_correction/{sample}.txt"
    log:
        "{sample}/logs/{sample}_error_correction.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    resources:
        mem = config.get("java_mem", JAVA_MEM)
    params:
        inputs=lambda wc, input: "in1={0},{2} in2={1}".format(*input) if paired_end else "in={0}".format(*input),
        outputs=lambda wc, output: "out1={0},{2} out2={1}".format(*output) if paired_end else "out={0}".format(*output)
    threads:
        config.get("threads", 1)
    shell:
        """
            {SHPFXM} tadpole.sh -Xmx{resources.mem}G \
               prealloc=1 \
               {params.inputs} \
               {params.outputs} \
               mode=correct \
               threads={threads} \
               ecc=t ecco=t 2>> {log}
        """


rule merge_pairs:
    input:
        expand("{{sample}}/assembly/reads/{{previous_steps}}_{fraction}.fastq.gz",
            fraction=multifile_fractions)
    output:
        temp(expand("{{sample}}/assembly/reads/{{previous_steps}}.merged_{fraction}.fastq.gz",
            fraction=multifile_fractions))
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", JAVA_MEM)
    conda:
        "%s/required_packages.yaml" % CONDAENV
    log:
        "{sample}/logs/{sample}_merge_pairs.log"
    benchmark:
        "logs/benchmarks/merge_pairs/{sample}.txt"
    shadow:
        "shallow"
    params:
        kmer = config.get("merging_k", MERGING_K),
        extend2 = config.get("merging_extend2", MERGING_EXTEND2),
        flags = config.get("merging_flags", MERGING_FLAGS)
    shell:
        """bbmerge.sh -Xmx{resources.mem}G threads={threads} \
               in1={input[0]} in2={input[1]} \
               outmerged={wildcards.sample}_merged_pairs.fastq.gz \
               outu={output[0]} outu2={output[1]} \
               {params.flags} k={params.kmer} \
               extend2={params.extend2} 2> {log}

           cat {wildcards.sample}_merged_pairs.fastq.gz {input[2]} \
               > {output[2]} 2>> {log}
        """




if config.get("assembler", "megahit") == "megahit":
    rule run_megahit:
        input:
            expand("{{sample}}/assembly/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
            fraction=multifile_fractions,assembly_preprocessing_steps=assembly_preprocessing_steps)
        output:
            temp("{sample}/assembly/{sample}_prefilter.contigs.fa")
        benchmark:
            "logs/benchmarks/assembly/{sample}.txt"
        shadow:
            "full"
        log:
             "{sample}/logs/{sample}_megahit.log"
        params:
            min_count = config.get("megahit_min_count", MEGAHIT_MIN_COUNT),
            k_min = config.get("megahit_k_min", MEGAHIT_K_MIN),
            k_max = config.get("megahit_k_max", MEGAHIT_K_MAX),
            k_step = config.get("megahit_k_step", MEGAHIT_K_STEP),
            merge_level = config.get("megahit_merge_level", MEGAHIT_MERGE_LEVEL),
            prune_level = config.get("megahit_prune_level", MEGAHIT_PRUNE_LEVEL),
            low_local_ratio = config.get("megahit_low_local_ratio", MEGAHIT_LOW_LOCAL_RATIO),
            min_contig_len = config.get("minimum_contig_length", MINIMUM_CONTIG_LENGTH),
            outdir = lambda wc, output: os.path.dirname(output[0]),
            inputs=lambda wc,input: "-1 {0} -2 {1} --read {2}".format(*input) if len(input)==2 else "--read {0}".format(*input)
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads:
            config.get("assembly_threads", ASSEMBLY_THREADS)
        resources:
            mem=config.get("assembly_memory", ASSEMBLY_MEMORY) #in GB
        shell:
            """{SHPFXM} megahit --continue \
                    {params.inputs} \
                   --tmp-dir {TMPDIR} \
                   --num-cpu-threads {threads} \
                   --k-min {params.k_min} \
                   --k-max {params.k_max} \
                   --k-step {params.k_step} \
                   --out-dir {params.outdir} \
                   --out-prefix {wildcards.sample}_prefilter \
                   --min-contig-len {params.min_contig_len} \
                   --min-count {params.min_count} \
                   --merge-level {params.merge_level} \
                   --prune-level {params.prune_level} \
                   --low-local-ratio {params.low_local_ratio} \
                   --memory {resources.mem}000000000  2> >(tee {log})
            """


    rule rename_megahit_output:
        input:
            "{sample}/assembly/{sample}_prefilter.contigs.fa"
        output:
            temp("{sample}/assembly/{sample}_raw_contigs.fasta")
        shell:
            "{SHPFXS} cp {input} {output}"

else:
    rule run_spades:
        input:
            expand("{{sample}}/assembly/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
                   fraction=multifile_fractions,
                   assembly_preprocessing_steps=assembly_preprocessing_steps)
        output:
            temp("{sample}/assembly/contigs.fasta")
        benchmark:
            "logs/benchmarks/assembly/{sample}.txt"
        params:
            inputs=lambda wc,input: "-1 {0} -2 {1} -s {2}".format(*input) if len(input) == 3 else "-s {0}".format(*input),
            k = config.get("spades_k", SPADES_K),
            outdir = lambda wc: "{sample}/assembly".format(sample=wc.sample),
            # FIXME: this does nothing
            error_correction= "" if False else "--only-assembler"
        log:
            "{sample}/logs/{sample}_spades.log"
        shadow:
            "full"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads:
            config.get("assembly_threads", ASSEMBLY_THREADS)
        resources:
            mem=config.get("assembly_memory", ASSEMBLY_MEMORY) #in GB
        shell:
            """
            spades.py --threads {threads} --memory {resources.mem} -o {params.outdir} --meta {params.inputs} \
            {params.error_correction} 2> >(tee {log}) """


    rule rename_spades_output:
        input:
            "{sample}/assembly/contigs.fasta"
        output:
            temp("{sample}/assembly/{sample}_raw_contigs.fasta")
        shell:
            "{SHPFXS} cp {input} {output}"


rule rename_contigs:
    # standardizes header labels within contig FASTAs
    input:
        "{sample}/assembly/{sample}_raw_contigs.fasta"
    output:
        "{sample}/assembly/{sample}_prefilter_contigs.fasta"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    shell:
        """rename.sh in={input} out={output} ow=t prefix={wildcards.sample}"""


rule calculate_prefiltered_contigs_stats:
    input:
        "{sample}/assembly/{sample}_prefilter_contigs.fasta"
    output:
        "{sample}/assembly/contig_stats/prefilter_contig_stats.txt"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        1
    resources:
        mem = config.get("java_mem", JAVA_MEM)
    shell:
        "{SHPFXS} stats.sh in={input} format=3 -Xmx{resources.mem}G > {output}"





rule calculate_prefiltered_contig_coverage_stats:
    input:
        unpack(get_quality_controlled_reads),
        fasta = "{sample}/assembly/{sample}_prefilter_contigs.fasta"
    output: # bbwrap gives output statistics only for single ended
       # bhist = "{sample}/assembly/contig_stats/prefilter_base_composition.txt",
       # bqhist = "{sample}/assembly/contig_stats/prefilter_box_quality.txt",
       # mhist = "{sample}/assembly/contig_stats/prefilter_mutation_rates.txt",
       # statsfile = "{sample}/assembly/contig_stats/prefilter_mapping_stats.txt",
        covstats = "{sample}/assembly/contig_stats/prefilter_coverage_stats.txt",
        sam= temp("{sample}/sequence_alignment/alignement_to_prefilter_contigs.sam.gz")
    benchmark:
        "logs/benchmarks/calculate_prefiltered_contig_coverage_stats/{sample}.txt"
    params:
        input= lambda wc,input : input_params_for_bbwrap(wc,input),
        interleaved = "auto" #lambda wc: "t" if config["samples"][wc.sample].get("paired", True) else "auto"
    log:
        "{sample}/assembly/logs/prefiltered_contig_coverage_stats.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", JAVA_MEM)
    shell:
        """{SHPFXM} bbwrap.sh nodisk=t ref={input.fasta} {params.input} fast=t \
               interleaved={params.interleaved} threads={threads} \
            -Xmx{resources.mem}G append out={output.sam} 2> {log}

            {SHPFXM} pileup.sh ref={input.fasta} in={output.sam} threads={threads} \
            -Xmx{resources.mem}G covstats={output.covstats} physcov 2>> {log}

        """


rule filter_by_coverage:
    input:
        fasta = "{sample}/assembly/{sample}_prefilter_contigs.fasta",
        covstats = "{sample}/assembly/contig_stats/prefilter_coverage_stats.txt"
    output:
        fasta = "{sample}/{sample}_contigs.fasta",
        removed_names = "{sample}/assembly/{sample}_discarded_contigs.fasta"
    params:
        minc = config.get("minimum_average_coverage", MINIMUM_AVERAGE_COVERAGE),
        minp = config.get("minimum_percent_covered_bases", MINIMUM_PERCENT_COVERED_BASES),
        minr = config.get("minimum_mapped_reads", MINIMUM_MAPPED_READS),
        minl = config.get("minimum_contig_length", MINIMUM_CONTIG_LENGTH),
        trim = config.get("contig_trim_bp", CONTIG_TRIM_BP)
    log:
        "{sample}/assembly/logs/filter_by_coverage.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        1
    resources:
        mem = config.get("java_mem", JAVA_MEM)
    shell:
        """{SHPFXS} filterbycoverage.sh in={input.fasta} \
               cov={input.covstats} \
               out={output.fasta} \
               outd={output.removed_names} \
               minc={params.minc} \
               minp={params.minp} \
               minr={params.minr} \
               minl={params.minl} \
               trim={params.trim} \
               -Xmx{resources.mem}G 2> {log}"""


rule align_reads_to_filtered_contigs:
    input:
        unpack(get_quality_controlled_reads),
        fasta = "{sample}/{sample}_contigs.fasta"
    output:
        # TODO: this should likely be temp
        sam = "{sample}/sequence_alignment/{sample}.sam",
        #bhist = "{sample}/assembly/contig_stats/postfilter_base_composition.txt",
        #bqhist = "{sample}/assembly/contig_stats/postfilter_box_quality.txt",
        #mhist = "{sample}/assembly/contig_stats/postfilter_mutation_rates.txt",
        #gchist = "{sample}/assembly/contig_stats/postfilter_gc_rates.txt",
        #statsfile = "{sample}/assembly/contig_stats/postfilter_mapping_stats.txt",
        basecov="{sample}/assembly/contig_stats/postfilter_base_coverage.txt.gz",
        covhist= "{sample}/assembly/contig_stats/postfilter_coverage_histogram.txt",
        covstats = "{sample}/assembly/contig_stats/postfilter_coverage_stats.txt",
        #unmapped=expand("{{sample}}/sequence_alignment/{{sample}}_unmapped_{fraction}.fastq.gz",fraction=interleaved_fractions)
    benchmark:
        "logs/benchmarks/align_reads_to_filtered_contigs/{sample}.txt"
    params:
        input= lambda wc,input : input_params_for_bbwrap(wc,input),
        #unmapped= lambda wc,output: ",".join(output.unmapped),
        interleaved = "auto", #lambda wc: "t" if config["samples"][wc.sample].get("paired", True) else "auto",
        maxsites = config.get("maximum_counted_map_sites", MAXIMUM_COUNTED_MAP_SITES)
    log:
        "{sample}/assembly/logs/contig_coverage_stats.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", JAVA_MEM)
    shell:
        """{SHPFXM} bbwrap.sh nodisk=t \
               ref={input.fasta} \
               {params.input} \
               trimreaddescriptions=t \
               out={output.sam} \
               mappedonly=t \
               threads={threads} \
               mdtag=t \
               xstag=fs \
               nmtag=t \
               sam=1.3 \
               local=t \
               ambiguous=all \
               interleaved={params.interleaved} \
               secondary=t \
               ssao=t \
               maxsites={params.maxsites} \
               -Xmx{resources.mem}G \
               append \
               2> {log}


            {SHPFXM} pileup.sh ref={input.fasta} in={output.sam} threads={threads} \
            -Xmx{resources.mem}G covstats={output.covstats} \
            hist={output.covhist} basecov={output.basecov} physcov 2>> {log}

            #samtools view -u -f4 {output.sam} | samtools bam2fq -s unmapped.se.fq - > unmapped.pe.fq

               """



if config.get("perform_genome_binning", True):
    rule make_maxbin_abundance_file:
        input:
            "{sample}/assembly/contig_stats/postfilter_coverage_stats.txt"
        output:
            "{sample}/genomic_bins/{sample}_contig_coverage.tsv"
        run:
            bb_cov_stats_to_maxbin(input[0], output[0])


    rule run_maxbin:
        input:
            fasta = "{sample}/{sample}_contigs.fasta",
            abundance = "{sample}/genomic_bins/{sample}_contig_coverage.tsv"
        output:
            # fastas will need to be dynamic if we do something with them at a later time
            summary = "{sample}/genomic_bins/{sample}.summary",
            marker = "{sample}/genomic_bins/{sample}.marker"
        benchmark:
            "logs/benchmarks/maxbin2/{sample}.txt"
        params:
            mi = config.get("maxbin_max_iteration", MAXBIN_MAX_ITERATION),
            mcl = config.get("maxbin_min_contig_length", MAXBIN_MIN_CONTIG_LENGTH),
            pt = config.get("maxbin_prob_threshold", MAXBIN_PROB_THRESHOLD),
            outdir = lambda wildcards, output: os.path.join(os.path.dirname(output.summary), wildcards.sample)
        log:
            "{sample}/logs/maxbin2.log"
        conda:
            "%s/optional_genome_binning.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        shell:
            """{SHPFXM} run_MaxBin.pl -contig {input.fasta} \
                   -abund {input.abundance} \
                   -out {params.outdir} \
                   -min_contig_length {params.mcl} \
                   -thread {threads} \
                   -prob_threshold {params.pt} \
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
            touched_output = "logs/checkm_init.txt",
            # init_checkm = "%s/hmms/checkm.hmm" % CHECKMDIR,
            bins = "{sample}/genomic_bins/{sample}.marker"
        output:
            "{sample}/genomic_bins/checkm/completeness.tsv"
        params:
            bin_dir = lambda wc, input: os.path.dirname(input.bins),
            output_dir = lambda wc, output: os.path.dirname(output[0])
        conda:
            "%s/optional_genome_binning.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        shell:
            """rm -r {params.output_dir} && \
               {SHPFXM} checkm lineage_wf \
                   --file {params.output_dir}/completeness.tsv \
                   --tab_table \
                   --quiet \
                   --extension fasta \
                   --threads {threads} \
                   {params.bin_dir} \
                   {params.output_dir}"""


    rule run_checkm_tree_qa:
        input:
            "{sample}/genomic_bins/checkm/completeness.tsv"
        output:
            "{sample}/genomic_bins/checkm/taxonomy.tsv"
        params:
            output_dir = "{sample}/genomic_bins/checkm"
        conda:
            "%s/optional_genome_binning.yaml" % CONDAENV
        shell:
            """{SHPFXS} checkm tree_qa \
                   --tab_table \
                   --out_format 2 \
                   --file {params.output_dir}/taxonomy.tsv \
                   {params.output_dir}"""


rule convert_sam_to_bam:
    input:
        "{file}.sam"
    output:
        "{file}.bam"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} samtools view \
               -@ {threads} \
               -bSh1 {input} | samtools sort \
                                   -m 1536M \
                                   -@ {threads} \
                                   -T {TMPDIR}/{wildcards.sample}_tmp \
                                   -o {output} \
                                   -O bam -"""


rule create_bam_index:
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        1
    shell:
        "{SHPFXS} samtools index {input}"


rule calculate_final_contigs_stats:
    input:
        "{sample}/{sample}_contigs.fasta"
    output:
        "{sample}/assembly/contig_stats/final_contig_stats.txt"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        1
    shell:
        "{SHPFXS} stats.sh in={input} format=3 > {output}"


rule run_prokka_annotation:
    input:
        "{sample}/{sample}_contigs.fasta"
    output:
        discrepancy = "{sample}/annotation/prokka/{sample}.err",
        faa = "{sample}/annotation/prokka/{sample}.faa",
        ffn = "{sample}/annotation/prokka/{sample}.ffn",
        fna = "{sample}/annotation/prokka/{sample}.fna",
        fsa = "{sample}/annotation/prokka/{sample}.fsa",
        gbk = "{sample}/annotation/prokka/{sample}.gbk",
        gff = "{sample}/annotation/prokka/{sample}.gff",
        log = "{sample}/annotation/prokka/{sample}.log",
        sqn = "{sample}/annotation/prokka/{sample}.sqn",
        tbl = "{sample}/annotation/prokka/{sample}.tbl",
        tsv = "{sample}/annotation/prokka/{sample}.tsv",
        txt = "{sample}/annotation/prokka/{sample}.txt"
    benchmark:
        "logs/benchmarks/prokka/{sample}.txt"
    params:
        outdir = lambda wc, output: os.path.dirname(output.faa),
        kingdom = config.get("prokka_kingdom", PROKKA_KINGDOM)
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} prokka --outdir {params.outdir} \
               --force \
               --prefix {wildcards.sample} \
               --locustag {wildcards.sample} \
               --kingdom {params.kingdom} \
               --metagenome \
               --cpus {threads} \
               {input}"""


rule update_prokka_tsv:
    input:
        "{sample}/annotation/prokka/{sample}.gff"
    output:
        "{sample}/annotation/prokka/{sample}_plus.tsv"
    shell:
        """{SHPFXS} atlas gff2tsv {input} {output}"""


rule convert_gff_to_gtf:
    input:
        "{sample}/annotation/prokka/{sample}.gff"
    output:
        "{sample}/annotation/prokka/{sample}.gtf"
    run:
        gff_to_gtf(input[0], output[0])


rule find_counts_per_region:
    input:
        gtf = "{sample}/annotation/prokka/{sample}.gtf",
        bam = "{sample}/sequence_alignment/{sample}.bam"
    output:
        summary = "{sample}/annotation/feature_counts/{sample}_counts.txt.summary",
        counts = "{sample}/annotation/feature_counts/{sample}_counts.txt"
    params:
        min_read_overlap = config.get("minimum_region_overlap", MINIMUM_REGION_OVERLAP),
        paired_mode = lambda wc: "-p" if config["samples"][wc.sample].get("paired", True) else "",
        multi_mapping = "-M" if config.get("count_multi_mapped_reads") else "",
        primary_only = "--primary" if config.get("primary_only", False) else ""
    log:
        "{sample}/logs/counts_per_region.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} featureCounts {params.paired_mode} \
               -F gtf \
               -T {threads} \
               {params.multi_mapping} \
               -t CDS \
               -g ID \
               -a {input.gtf} \
               -o {output.counts} \
               {input.bam} 2> {log}"""


rule run_diamond_blastp:
    input:
        fasta = "{sample}/annotation/prokka/{sample}.faa",
        db = config["diamond_db"]
    output:
        "{sample}/annotation/refseq/{sample}_hits.tsv"
    benchmark:
        "logs/benchmarks/run_diamond_blastp/{sample}.txt"
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
        """{SHPFXM} diamond blastp \
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


rule add_contig_metadata:
    input:
        hits = "{sample}/annotation/refseq/{sample}_hits.tsv",
        gff = "{sample}/annotation/prokka/{sample}.gff"
    output:
        temp("{sample}/annotation/refseq/{sample}_hits_plus.tsv")
    shell:
        "{SHPFXS} atlas munge-blast {input.hits} {input.gff} {output}"


rule sort_munged_blast_hits:
    # ensure blast hits are grouped by contig, ORF, and then decreasing by bitscore
    input:
        "{sample}/annotation/refseq/{sample}_hits_plus.tsv"
    output:
        "{sample}/annotation/refseq/{sample}_hits_plus_sorted.tsv"
    shell:
        "{SHPFXS} sort -k1,1 -k2,2 -k13,13rn {input} > {output}"


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
        """{SHPFXS} atlas refseq \
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
               {output}"""


if config.get("perform_genome_binning", True):
    rule merge_sample_tables:
        input:
            prokka = "{sample}/annotation/prokka/{sample}_plus.tsv",
            refseq = "{sample}/annotation/refseq/{sample}_tax_assignments.tsv",
            counts = "{sample}/annotation/feature_counts/{sample}_counts.txt",
            completeness = "{sample}/genomic_bins/checkm/completeness.tsv",
            taxonomy = "{sample}/genomic_bins/checkm/taxonomy.tsv"
        output:
            "{sample}/{sample}_annotations.txt"
        params:
            fastas = lambda wc: " --fasta ".join(glob("{sample}/genomic_bins/{sample}.*.fasta".format(sample=wc.sample)))
        shell:
            "{SHPFXS} atlas merge-tables \
                 --counts {input.counts} \
                 --completeness {input.completeness} \
                 --taxonomy {input.taxonomy} \
                 --fasta {params.fastas} \
                 {input.prokka} \
                 {input.refseq} \
                 {output}"


else:
    rule merge_sample_tables:
        input:
            prokka = "{sample}/annotation/prokka/{sample}_plus.tsv",
            refseq = "{sample}/annotation/refseq/{sample}_tax_assignments.tsv",
            counts = "{sample}/annotation/feature_counts/{sample}_counts.txt"
        output:
            "{sample}/{sample}_annotations.txt"
        shell:
            "{SHPFXS} atlas merge-tables \
                 --counts {input.counts} \
                 {input.prokka} \
                 {input.refseq} \
                 {output}"


# rule assembly_report:
#
#     input:
#         contig_stats = "{sample}/assembly/contig_stats/final_contig_stats.txt",
#         base_comp = "{sample}/assembly/contig_stats/postfilter_base_composition.txt"
#         # css = os.path.join(workflow.basedir, "resources", "report.css")
#     output:
#         html = "{sample}/{sample}_assembly_README.html"
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
