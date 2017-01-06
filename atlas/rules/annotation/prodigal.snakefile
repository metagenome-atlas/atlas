import re


def gff_to_gtf(gff_in, gtf_out):
    t = re.compile(r'ID=[0-9]+_([0-9]+);')
    with open(gtf_out, "w") as fh, open(gff_in) as gff:
        for line in gff:
            if line.startswith("#"): continue
            toks = line.strip().split("\t")
            orf = t.findall(toks[-1])[0]
            gene_id = toks[0] + "_" + orf
            toks[-1] = 'gene_id "%s"; %s' % (gene_id, toks[-1])
            print(*toks, sep="\t", file=fh)


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
