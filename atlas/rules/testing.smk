
from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
NCBI = NCBIRemoteProvider(email="someone@example.com") # email required by NCBI to prevent abuse


Genomes={
    "Mycoplasma":"NC_019552.1.fasta",#hyorhinis
    "Streptococcus":"NC_017581.1.fasta", #thermophilus
    #"Ureaplasma": "NC_011374.1.fasta" #urealyticum
}


rule all:
    input:
        expand("{genome}_{fraction}.fastq.gz",
               genome=Genomes.keys(),
               fraction=['R1','R2'])

# rule download_genome:
#     input:
#         lambda wildcards: NCBI.remote(Genomes[wildcards.genome], db="nuccore")
#     output:
#         temp("{genome}.fasta")
#     shell:
#         "cp {input} {output}"
#
# rule get_metagenome:
#     input:
#         expand(rules.download_genome.output,genome=Genomes.keys())
#     output:
#         temp('Metagenome.fasta')
#     shell:
#         "cat {input} > {output}"

rule gen_random_reads:
    input:
        lambda wildcards: NCBI.remote(Genomes[wildcards.genome], db="nuccore")
    output:
        temp("{genome}.fastq.gz") # interleaved
    params:
        reads= config.get('reads',5e4),
        maxlength=120,
        metagenome=lambda wc: 't' if wc.genome=='Metagenome' else 'f'
    shell:
        " randomreads.sh "
        " ref={input} "
        " paired=t illuminanames=t "
        " metagenome={params.metagenome} "
        " maxlength={params.maxlength} reads={params.reads}"
        " out={output}"

rule de_interleave:
    input:
        rules.gen_random_reads.output
    output:
        expand("{{genome}}_{fraction}.fastq.gz",
               fraction=['R1','R2'])
    shell:
        "reformat.sh in={input} out={output[0]} out2={output[1]}"
