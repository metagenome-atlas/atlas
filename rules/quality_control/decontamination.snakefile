def get_handle_ribosomal_rna_input(wildcards):
    inputs = []
    data_type = config["samples"][wildcards.sample].get("type", "metagenome").lower()
    if data_type == "metagenome":
        return ["{sample}/quality_control/decontamination/{sample}_clean.fastq.gz".format(sample=wildcards.sample),
                "{sample}/quality_control/decontamination/{sample}_rRNA.fastq.gz".format(sample=wildcards.sample)]
    else:
        return ["{sample}/quality_control/decontamination/{sample}_clean.fastq.gz".format(sample=wildcards.sample)]


rule decontamination:
    input:
        "{sample}/quality_control/error_correction/{sample}_pe.fastq.gz"
    output:
        dbs = ["{sample}/quality_control/decontamination/{sample}_%s.fastq.gz" % db for db in list(config["preprocessing"]["contamination"]["references"].keys())],
        stats = "{sample}/quality_control/decontamination/{sample}_refstats.txt",
        clean = "{sample}/quality_control/decontamination/{sample}_clean.fastq.gz"
    params:
        refs_in = " ".join(["ref_%s=%s" % (n, fa) for n, fa in config["preprocessing"]["contamination"]["references"].items()]),
        refs_out = lambda wc: " ".join(["out_{ref}={sample}/quality_control/decontamination/{sample}_{ref}.fastq.gz".format(ref=n, sample=wc.sample) for n in list(config["preprocessing"]["contamination"]["references"].keys())]),
        maxindel = config["preprocessing"]["contamination"].get("maxindel", 20),
        minratio = config["preprocessing"]["contamination"].get("minratio", 0.65),
        minhits = config["preprocessing"]["contamination"].get("minhits", 1),
        ambiguous = config["preprocessing"]["contamination"].get("ambiguous", "best"),
        k = config["preprocessing"]["contamination"].get("k", 13)
    log:
        "{sample}/logs/{sample}_decontamination.log"
    threads:
        config.get("threads", 1)
    shell:
        """{SHPFXM} bbsplit.sh {params.refs_in} in={input} outu={output.clean} \
               {params.refs_out} maxindel={params.maxindel} minratio={params.minratio} \
               minhits={params.minhits} ambiguous={params.ambiguous} refstats={output.stats}\
               threads={threads} k={params.k} local=t 2> {log}"""


rule handle_ribosomal_rna:
    input:
        get_handle_ribosomal_rna_input
    output:
        "{sample}/quality_control/decontamination/{sample}_pe.fastq.gz"
    threads:
        1
    shell:
        "{SHPFXS} cat {input} > {output}"
