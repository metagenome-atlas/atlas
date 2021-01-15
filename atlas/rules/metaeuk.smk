
# this is the target rule
# To acitvte this workflow run 'atlas run None eukaryotes'
rule eukaryotes:
    input:
        expand("{sample}/taxonomy/mmseqs_output",sample=SAMPLES)

MMseqs_DB_folder=os.path.join(DBDIR,'MMSeqs')


rule download_mmseqs:
    output:
        f"{MMseqs_DB_folder}/swissprot"
    conda:
        "../envs/metaeuk.yaml"
    shell:
        'echo "Not yet implemented!"'


rule contig_taxonomy:
    input:
        contigs="{sample}/{sample}_contigs.fasta",
        db= f"{MMseqs_DB_folder}/swissprot"
    output:
        directory("{sample}/taxonomy/mmseqs_output"),
    params:
        tmpdir= os.path.join(config['tmpdir'],"mmseqs"),
    log:
        "{sample}/logs/taxonomy/mmseqs_taxonomy.log"
    threads:
        config["threads"]
    resources:
        mem= config['mem'], #or large_mem
        time= config['time']
    conda:
        "../envs/metaeuk.yaml"
    shell:
        "mmseqs easy-taxonomy "
        " {input.contigs} {input.db}"
        " {output} {params.tmpdir} "
        " --search-type 2"
