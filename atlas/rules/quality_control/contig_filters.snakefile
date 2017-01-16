rule dirty_contigs_stats:
    input:
        "{sample}/{ASSEMBLER}/{sample}_prefilter_contigs.fasta"
    output:
        "{sample}/{ASSEMBLER}/stats/prefilter_contig_stats.txt"
    threads:
        1
    shell:
        "{SHPFXS} stats.sh in={input} format=3 > {output}"


rule dirty_contig_coverage_stats:
    input:
        fasta = "{sample}/{ASSEMBLER}/{sample}_prefilter_contigs.fasta",
        fastq = "{sample}/quality_control/decontamination/{sample}_pe.fastq.gz"
    output:
        bhist = "{sample}/{ASSEMBLER}/stats/prefilter_base_composition.txt",
        bqhist = "{sample}/{ASSEMBLER}/stats/prefilter_box_quality.txt",
        mhist = "{sample}/{ASSEMBLER}/stats/prefilter_mutation_rates.txt",
        statsfile = "{sample}/{ASSEMBLER}/stats/prefilter_mapping_stats.txt",
        covstats = "{sample}/{ASSEMBLER}/stats/prefilter_coverage_stats.txt"
    log:
        "{sample}/{ASSEMBLER}/logs/dirty_contig_coverage_stats.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbmap.sh nodisk=t ref={input.fasta} in={input.fastq} fast=t \
               threads={threads} bhist={output.bhist} bqhist={output.bqhist} mhist={output.mhist} \
               statsfile={output.statsfile} covstats={output.covstats} 2> {log}"""


rule filter_by_coverage:
    input:
        fasta = "{sample}/{ASSEMBLER}/{sample}_prefilter_contigs.fasta",
        covstats = "{sample}/{ASSEMBLER}/stats/prefilter_coverage_stats.txt"
    output:
        fasta = "{sample}/{ASSEMBLER}/{sample}_contigs.fasta",
        removed_names = "{sample}/{ASSEMBLER}/{sample}_discarded_contigs.txt"
    params:
        minc = config["assembly"].get("minc", 5),
        minp = config["assembly"].get("minp", 40),
        minr = config["assembly"].get("minr", 0),
        minl = config["assembly"].get("minl", 1),
        trim = config["assembly"].get("trim", 0)
    log:
        "{sample}/{ASSEMBLER}/logs/filter_by_coverage.log"
    threads:
        1
    shell:
        """{SHPFXS} filterbycoverage.sh in={input.fasta} cov={input.covstats} out={output.fasta} \
               outd={output.removed_names} minc={params.minc} minp={params.minp} \
               minr={params.minr} minl={params.minl} trim={params.trim} 2> {log}"""


rule contig_coverage_stats:
    input:
        fasta = "{sample}/{ASSEMBLER}/{sample}_contigs.fasta",
        fastq = "{sample}/quality_control/decontamination/{sample}_pe.fastq.gz"
    output:
        sam = temp("{sample}/{ASSEMBLER}/annotation/{sample}.sam"),
        bhist = "{sample}/{ASSEMBLER}/stats/postfilter_base_composition.txt",
        bqhist = "{sample}/{ASSEMBLER}/stats/postfilter_box_quality.txt",
        mhist = "{sample}/{ASSEMBLER}/stats/postfilter_mutation_rates.txt",
        gchist = "{sample}/{ASSEMBLER}/stats/postfilter_gc_rates.txt",
        statsfile = "{sample}/{ASSEMBLER}/stats/postfilter_mapping_stats.txt",
        covstats = "{sample}/{ASSEMBLER}/stats/postfilter_coverage_stats.txt"
    log:
        "{sample}/{ASSEMBLER}/logs/contig_coverage_stats.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbmap.sh nodisk=t ref={input.fasta} in={input.fastq} trimreaddescriptions=t \
               out={wildcards.sample}/{ASSEMBLER}/annotation/{wildcards.sample}.sam \
               mappedonly=t threads={threads} bhist={output.bhist} bqhist={output.bqhist} \
               mhist={output.mhist} gchist={output.gchist} statsfile={output.statsfile} \
               covstats={output.covstats} mdtag=t xstag=fs nmtag=t sam=1.3 2> {log}"""


rule sam_to_bam:
   input:
       "{sample}/{ASSEMBLER}/annotation/{sample}.sam"
   output:
       "{sample}/{ASSEMBLER}/annotation/{sample}.bam"
   threads:
       config.get("threads", 1)
   shell:
       """{SHPFXM} samtools view -@ {threads} -bSh1 {input} | samtools sort -@ {threads} -T {TMPDIR}/{wildcards.sample}_tmp -o {output} -O bam -"""


rule create_bam_index:
   input:
       "{sample}/{ASSEMBLER}/annotation/{sample}.bam"
   output:
       "{sample}/{ASSEMBLER}/annotation/{sample}.bam.bai"
   threads:
       1
   shell:
       "{SHPFXS} samtools index {input}"


rule final_contigs_stats:
   input:
       "{sample}/{ASSEMBLER}/{sample}_contigs.fasta"
   output:
       "{sample}/{ASSEMBLER}/stats/final_contig_stats.txt"
   threads:
       1
   shell:
       "{SHPFXS} stats.sh in={input} format=3 > {output}"
