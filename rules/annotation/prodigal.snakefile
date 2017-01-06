from atlas.utils import gff_to_gtf


rule prodigal_orfs:
    input:
        "{sample}/%s/{sample}_contigs.fasta" % ASSEMBLER
    output:
        prot = "{sample}/annotation/orfs/{sample}.faa",
        nuc = "{sample}/annotation/orfs/{sample}.fna",
        gff = "{sample}/annotation/orfs/{sample}.gff"
    params:
        g = config["annotation"].get("translation_table", "11")
    threads:
        1
    shell:
        """{SHPFXS} prodigal -i {input} -o {output.gff} -f gff -a {output.prot} -d {output.nuc} \
               -g {params.g} -p meta"""


rule gff_to_gtf:
    input:
        "{sample}/annotation/orfs/{sample}.gff"
    output:
        "{sample}/annotation/orfs/{sample}.gtf"
    run:
        gff_to_gtf(input[0], output[0])
