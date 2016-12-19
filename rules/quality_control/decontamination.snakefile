def get_handle_ribosomal_rna_input(wildcards):
    inputs = []
    data_type = config["samples"][wildcards.sample].get("type", "metagenome").lower()
    if data_type == "metagenome":
        return ["results/{eid}/{sample}/quality_control/decontamination/{sample}_clean.fastq.gz".format(eid=wildcards.eid, sample=wildcards.sample),
                "results/{eid}/{sample}/quality_control/decontamination/{sample}_rRNA.fastq.gz".format(eid=wildcards.eid, sample=wildcards.sample)]
    else:
        return ["results/{eid}/{sample}/quality_control/decontamination/{sample}_clean.fastq.gz".format(eid=wildcards.eid, sample=wildcards.sample)]


rule decontamination:
    input:
        "results/{eid}/{sample}/quality_control/error_correction/{sample}_pe.fastq.gz"
    output:
        dbs = ["results/{eid}/{sample}/quality_control/decontamination/{sample}_%s.fastq.gz" % db for db in list(config["preprocessing"]["contamination"]["references"].keys())],
        stats = "results/{eid}/{sample}/quality_control/decontamination/{sample}_refstats.txt",
        clean = "results/{eid}/{sample}/quality_control/decontamination/{sample}_clean.fastq.gz"
    params:
        refs_in = " ".join(["ref_%s=%s" % (n, fa) for n, fa in config["preprocessing"]["contamination"]["references"].items()]),
        refs_out = lambda wc: " ".join(["out_%s=results/%s/%s/quality_control/decontamination/%s_%s.fastq.gz" % (n, wc.eid, wc.sample, wc.sample, n) for n in list(config["preprocessing"]["contamination"]["references"].keys())]),
        maxindel = config["preprocessing"]["contamination"].get("maxindel", 20),
        minratio = config["preprocessing"]["contamination"].get("minratio", 0.65),
        minhits = config["preprocessing"]["contamination"].get("minhits", 1),
        ambiguous = config["preprocessing"]["contamination"].get("ambiguous", "best"),
        k = config["preprocessing"]["contamination"].get("k", 15)
    log:
        "results/{eid}/{sample}/logs/{sample}_decontamination.log"
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
        "results/{eid}/{sample}/quality_control/decontamination/{sample}_pe.fastq.gz"
    threads:
        1
    shell:
        "{SHPFXS} cat {input} > {output}"
