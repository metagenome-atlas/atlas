
rule all:

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

rule coassembly_contig_coverage:

rule coassembly_contig_filter:

rule coassembly_postfilter_stats:

rule coassembly_sample_mapping:

rule coassembly_sam_to_bam:

rule coassembly_index_bam:

rule coassembly_prodigal:

gff to gtf
verse
diamond across databases
munging

    include: "rules/annotation/diamond.snakefile"
    include: "rules/annotation/prodigal.snakefile"
    include: "rules/annotation/verse.snakefile"
    include: "rules/annotation/munging.snakefile"


cat shewanella/quality_control/normalization_k19_t100/shewanella_pe.fastq.gz cytophaga/quality_control/normalization_k19_t100/cytophaga_pe.fastq.gz flavobacterium/quality_control/normalization_k19_t100/flavobacterium_pe.fastq.gz > coassemblies/test-coa/all_normalized_reads/test-coa_pe.fastq.gz
mkdir -p coassemblies/test-coa/all_decontamination_reads
cat shewanella/quality_control/decontamination/shewanella_pe.fastq.gz cytophaga/quality_control/decontamination/cytophaga_pe.fastq.gz flavobacterium/quality_control/decontamination/flavobacterium_pe.fastq.gz > coassemblies/test-coa/all_decontamination_reads/test-coa_pe.fastq.gz


bbnorm.sh in=coassemblies/test-coa/all_normalized_reads/test-coa_pe.fastq.gz out=coassemblies/test-coa/normalization_k19_t100/test-coa_pe.fastq.gz k=19 t=100 minkmers=8 prefilter=t threads=24
megahit --num-cpu-threads 24 --12 coassemblies/test-coa/normalization_k19_t100/test-coa_pe.fastq.gz --continue --k-min 21 --k-max 121 --k-step 20 --out-dir coassemblies/test-coa/megahit_21_121_20_normalization_k19_t100 --out-prefix test-coa_prefilter
mv coassemblies/test-coa/megahit_21_121_20_normalization_k19_t100/test-coa_prefilter.contigs.fa coassemblies/test-coa/megahit_21_121_20_normalization_k19_t100/test-coa_prefilter_contigs.fa
mkdir -p coassemblies/test-coa/megahit_21_121_20_normalization_k19_t100/stats/
stats.sh in=coassemblies/test-coa/megahit_21_121_20_normalization_k19_t100/test-coa_prefilter_contigs.fa format=3 > coassemblies/test-coa/megahit_21_121_20_normalization_k19_t100/stats/prefilter_contig_stats.txt
bbmap.sh nodisk=t ref=coassemblies/test-coa/megahit_21_121_20_normalization_k19_t100/test-coa_prefilter_contigs.fa in=coassemblies/test-coa/all_decontamination_reads/test-coa_pe.fastq.gz fast=t threads=24 covstats=coassemblies/test-coa/megahit_21_121_20_normalization_k19_t100/stats/prefilter_coverage_stats.txt

filterbycoverage.sh in=coassemblies/test-coa/megahit_21_121_20_normalization_k19_t100/test-coa_prefilter_contigs.fa cov=coassemblies/test-coa/megahit_21_121_20_normalization_k19_t100/stats/prefilter_coverage_stats.txt out=coassemblies/test-coa/megahit_21_121_20_normalization_k19_t100/test-coa_contigs.fa minc=5 minp=40 minr=0 minl=200 trim=50
stats.sh in=coassemblies/test-coa/megahit_21_121_20_normalization_k19_t100/test-coa_contigs.fa format=3 > coassemblies/test-coa/megahit_21_121_20_normalization_k19_t100/stats/final_contig_stats.txt
bbmap.sh nodisk=t ref=coassemblies/test-coa/megahit_21_121_20_normalization_k19_t100/test-coa_contigs.fa in=shewanella/quality_control/decontamination/shewanella_pe.fastq.gz threads=24 trimreaddescriptions=t out=coassemblies/test-coa/shewanella/annotation/shewanella.sam mappedonly=t mdtag=t xstag=fs nmtag=t sam=1.3
bbmap.sh nodisk=t ref=coassemblies/test-coa/megahit_21_121_20_normalization_k19_t100/test-coa_contigs.fa in=flavobacterium/quality_control/decontamination/flavobacterium_pe.fastq.gz threads=24 trimreaddescriptions=t out=coassemblies/test-coa/flavobacterium/annotation/flavobacterium.sam mappedonly=t mdtag=t xstag=fs nmtag=t sam=1.3
bbmap.sh nodisk=t ref=coassemblies/test-coa/megahit_21_121_20_normalization_k19_t100/test-coa_contigs.fa in=cytophaga/quality_control/decontamination/cytophaga_pe.fastq.gz threads=24 trimreaddescriptions=t out=coassemblies/test-coa/cytophaga/annotation/cytophaga.sam mappedonly=t mdtag=t xstag=fs nmtag=t sam=1.3
samtools view -@ 24 -bSh1 coassemblies/test-coa/shewanella/annotation/shewanella.sam | samtools sort -@ 24 -T shewanella_tmp -o coassemblies/test-coa/shewanella/annotation/shewanella.bam -O bam -
samtools view -@ 24 -bSh1 coassemblies/test-coa/flavobacterium/annotation/flavobacterium.sam | samtools sort -@ 24 -T flavobacterium_tmp -o coassemblies/test-coa/flavobacterium/annotation/flavobacterium.bam -O bam -
samtools view -@ 24 -bSh1 coassemblies/test-coa/cytophaga/annotation/cytophaga.sam | samtools sort -@ 24 -T cytophaga_tmp -o coassemblies/test-coa/cytophaga/annotation/cytophaga.bam -O bam -
samtools index coassemblies/test-coa/shewanella/annotation/shewanella.bam
samtools index coassemblies/test-coa/flavobacterium/annotation/flavobacterium.bam
samtools index coassemblies/test-coa/cytophaga/annotation/cytophaga.bam


mapping of each sample back to this coassembly
generate the counts per sample
compile these counts into a single table with metric, sample count
