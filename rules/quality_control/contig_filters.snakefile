rule dirty_contigs_stats:
    input:
        "{sample}/%s/{sample}_prefilter_contigs.fasta" % ASSEMBLER
    output:
        "{sample}/%s/stats/prefilter_contig_stats.txt" % ASSEMBLER
    threads:
        1
    shell:
        "{SHPFXS} stats.sh in={input} format=3 > {output}"


rule dirty_contig_coverage_stats:
    input:
        fasta = "{sample}/%s/{sample}_prefilter_contigs.fasta" % ASSEMBLER,
        fastq = "{sample}/quality_control/decontamination/{sample}_pe.fastq.gz"
    output:
        bhist = "{sample}/%s/stats/prefilter_base_composition.txt" % ASSEMBLER,
        bqhist = "{sample}/%s/stats/prefilter_box_quality.txt" % ASSEMBLER,
        mhist = "{sample}/%s/stats/prefilter_mutation_rates.txt" % ASSEMBLER,
        statsfile = "{sample}/%s/stats/prefilter_mapping_stats.txt" % ASSEMBLER,
        covstats = "{sample}/%s/stats/prefilter_coverage_stats.txt" % ASSEMBLER
    log:
        "{sample}/logs/dirty_contig_coverage_stats.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbmap.sh nodisk=t ref={input.fasta} in={input.fastq} fast=t \
               threads={threads} bhist={output.bhist} bqhist={output.bqhist} mhist={output.mhist} \
               statsfile={output.statsfile} covstats={output.covstats} 2> {log}"""


rule filter_by_coverage:
    input:
        fasta = "{sample}/%s/{sample}_prefilter_contigs.fasta" % ASSEMBLER,
        covstats = "{sample}/%s/stats/prefilter_coverage_stats.txt" % ASSEMBLER
    output:
        fasta = "{sample}/%s/{sample}_contigs.fasta" % ASSEMBLER,
        removed_names = "{sample}/%s/{sample}_discarded_contigs.txt" % ASSEMBLER
    params:
        minc = config["assembly"].get("minc", 5),
        minp = config["assembly"].get("minp", 40),
        minr = config["assembly"].get("minr", 0),
        minl = config["assembly"].get("minl", 1),
        trim = config["assembly"].get("trim", 0)
    log:
        "{sample}/logs/filter_by_coverage.log"
    threads:
        1
    shell:
        """{SHPFXS} filterbycoverage.sh in={input.fasta} cov={input.covstats} out={output.fasta} \
               outd={output.removed_names} minc={params.minc} minp={params.minp} \
               minr={params.minr} minl={params.minl} trim={params.trim} 2> {log}"""


rule contig_coverage_stats:
    input:
        fasta = "{sample}/%s/{sample}_contigs.fasta" % ASSEMBLER,
        fastq = "{sample}/quality_control/decontamination/{sample}_pe.fastq.gz"
    output:
        sam = temp("{sample}/annotation/{sample}.sam"),
        # bai = "{sample}/annotation/{sample}.bam.bai",
        bhist = "{sample}/%s/stats/postfilter_base_composition.txt" % ASSEMBLER,
        bqhist = "{sample}/%s/stats/postfilter_box_quality.txt" % ASSEMBLER,
        mhist = "{sample}/%s/stats/postfilter_mutation_rates.txt" % ASSEMBLER,
        gchist = "{sample}/%s/stats/postfilter_gc_rates.txt" % ASSEMBLER,
        statsfile = "{sample}/%s/stats/postfilter_mapping_stats.txt" % ASSEMBLER,
        covstats = "{sample}/%s/stats/postfilter_coverage_stats.txt" % ASSEMBLER
    log:
        "{sample}/logs/contig_coverage_stats.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbmap.sh nodisk=t ref={input.fasta} in={input.fastq} trimreaddescriptions=t \
               out={wildcards.sample}/annotation/{wildcards.sample}.sam \
               mappedonly=t threads={threads} bhist={output.bhist} bqhist={output.bqhist} \
               mhist={output.mhist} gchist={output.gchist} statsfile={output.statsfile} \
               covstats={output.covstats} mdtag=t xstag=fs nmtag=t sam=1.3 2> {log}"""


rule sam_to_bam:
   input:
       "{sample}/annotation/{sample}.sam"
   output:
       "{sample}/annotation/{sample}.bam"
   threads:
       config.get("threads", 1)
   shell:
       """{SHPFXM} samtools view -@ {threads} -bSh1 {input} | samtools sort -@ {threads} -T {TMPDIR}/{wildcards.sample}_tmp -o {output} -O bam -"""


rule create_bam_index:
   input:
       "{sample}/annotation/{sample}.bam"
   output:
       "{sample}/annotation/{sample}.bam.bai"
   shell:
       "{SHPFXS} samtools index {input}"


rule final_contigs_stats:
   input:
       "{sample}/%s/{sample}_contigs.fasta" % ASSEMBLER
   output:
       "{sample}/%s/stats/final_contig_stats.txt" % ASSEMBLER
   threads:
       1
   shell:
       "{SHPFXS} stats.sh in={input} format=3 > {output}"
