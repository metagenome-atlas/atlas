localrules:
    move_qc_reads,
    copy_assembly,


# Rules that are usefull temporarily to update to new version of atlas




rule copy_assembly:
    input:
        "{sample}/{sample}_contigs.fasta",
    output:
        "Assembly/fasta/{sample}.fasta",
    shell:
        "cp {input} {output}"
