SEMIBIN_DATA_PATH = os.path.join(DBDIR, "SemiBin_GTDB")


rule semibin_download_gtdb:
    output:
        directory(SEMIBIN_DATA_PATH),
    log:
        "log/download/Semibin.txt",
    conda:
        "../envs/semibin.yaml"
    threads: 1
    shell:
        "SemiBin --reference-db {output} 2> {log}"


rule semibin_predict_taxonomy:
    input:
        fasta=rules.combine_contigs.output,
        db=SEMIBIN_DATA_PATH,
    output:
        "Crossbinning/SemiBin/inconsitent_taxonomy.tsv",
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["default"],
    log:
        "log/semibin/predict_taxonomy.log",
    benchmark:
        "log/benchmarks/semibin/predict_taxonomy.tsv"
    params:
        output_dir=lambda wc, output: os.path.dirname(output[0]),
        name=lambda wc, output: os.path.basename(output[0]),
    shell:
        "SemiBin predict_taxonomy "
        " --input-fasta {input.fasta} "
        " --outdir {params.output_dir} "
        " --cannot-name {params.name} "
        " --threads {threads} "
        " --reference-db {input.db} "
        "2> {log}"


rule semibin_generate_data_multi:
    input:
        fasta=rules.combine_contigs.output,
        bams=expand(rules.sort_bam.output, sample=SAMPLES),
    output:
        "Crossbinning/SemiBin/data.csv",
        "Crossbinning/SemiBin/data_split.csv",
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["default"],
    log:
        "log/semibin/generate_data_multi.log",
    benchmark:
        "log/benchmarks/semibin/generate_data_multi.tsv"
    params:
        output_dir=lambda wc, output: os.path.dirname(output[0]),
        separator=":",
    shell:
        "SemiBin generate_data_multi "
        " --input-fasta {input.fasta} "
        " --input-bam {input.bams} "
        " --outdir {params.output_dir} "
        " --threads {threads} "
        " --separator {params.separator}"
        "2> {log}"


rule semibin_train:
    input:
        fasta=rules.combine_contigs.output,
        bams=expand(rules.sort_bam.output, sample=SAMPLES),
        data="Crossbinning/SemiBin/data.csv",
        data_split="Crossbinning/SemiBin/data_split.csv",
        cannot_link=rules.semibin_predict_taxonomy.output[0],
    output:
        "Crossbinning/SemiBin/model.h5",
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["default"],
    log:
        "log/semibin/train.log",
    benchmark:
        "log/benchmarks/semibin/train.tsv"
    params:
        output_dir=lambda wc, output: os.path.dirname(output[0]),
        extra=" --epoches 20",
    shell:
        "SemiBin train "
        " --input-fasta {input.fasta} "
        " --input-bam {input.bams} "
        " --outdir {params.output_dir} "
        " --threads {threads} "
        " --data {input.data} "
        " --data_split {input.data_split} "
        " --cannot-link {input.cannot_link} "
        " {extra} "
        "2> {log}"


rule semibin:
    input:
        fasta=rules.combine_contigs.output,
        bams=expand(rules.sort_bam.output, sample=SAMPLES),
        data="Crossbinning/SemiBin/data.csv",
        model=rules.semibin_train.output[0],
    output:
        touch("Crossbinning/SemiBin/finished_binning"),
    conda:
        "../envs/semibin.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        time=config["runtime"]["default"],
    log:
        "log/semibin/bin.log",
    benchmark:
        "log/benchmarks/semibin/bin.tsv"
    params:
        output_dir=lambda wc, output: os.path.dirname(output[0]),
        extra=" --minfasta-kbs 200 --recluster --max-node 1 --max-edges 200 ",
    shell:
        "SemiBin train "
        " --input-fasta {input.fasta} "
        " --input-bam {input.bams} "
        " --outdir {params.output_dir} "
        " --threads {threads} "
        " --data {input.data} "
        " --model {input.model} "
        " {params.extra} "
        "2> {log}"


# alternative to pretrained model --environment: Environment for the built-in model(human_gut/dog_gut/ocean).”
