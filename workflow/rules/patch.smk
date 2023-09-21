localrules: move_qc_reads, copy_assembly

# Rules that are usefull temporarily to update to new version of atlas


ruleorder: move_qc_reads > qcreads
rule move_qc_reads:
    input:
        "{sample}/sequence_quality_control/{sample}_clean_{fraction}.fastq.gz"
    output:
        "QC/reads/{sample}_{fraction}.fastq.gz"
    shell:
        "mv {input} {output}"

rule copy_assembly:
    input:
        "{sample}/{sample}_contigs.fasta"
    output:
        "Assembly/fasta/{sample}.fasta"
    shell:
        "cp {input} {output}"

