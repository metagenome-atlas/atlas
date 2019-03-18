


from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
NCBI = NCBIRemoteProvider(email="someone@example.com") # email required by NCBI to prevent abuse

rule all:
    input:
        "size.txt"

Genomes={
    "Mycoplasma_hyorhinis":"NC_019552.1.fasta"


}


rule download_genomes:
    input:
        NCBI.remote("KY785484.1.fasta", db="nuccore")
    output:
        "size.txt"
    run:
        shell("wc -c {input} > {output}")
