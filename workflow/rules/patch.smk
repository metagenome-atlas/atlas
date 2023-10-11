localrules:
    copy_assembly,


# Rules that are usefull temporarily to update to new version of atlas

ruleorder: copy_assembly > finalize_contigs

rule copy_assembly:
    input:
        "{sample}/{sample}_contigs.fasta",
    output:
        "Assembly/fasta/{sample}.fasta",
    shell:
        "cp {input} {output}"
