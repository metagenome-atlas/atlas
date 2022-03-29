wildcard_constraints:
    sra_run="[S,E,D]RR[0-9]+",


localrules:
    prefetch,


SRA_read_fractions = ["_1", "_2"] if PAIRED_END else [""]
SRA_SUBDIR_RUN = "SRA/Runs"


rule prefetch:
    output:
        sra=temp(touch(SRA_SUBDIR_RUN + "/{sra_run}/{sra_run}_downloaded")),
        # not givins sra file as output allows for continue from the same download
    params:
        outdir=SRA_SUBDIR_RUN,  # prefetch creates file in subfolder with run name automatically
    log:
        "logs/SRAdownload/prefetch/{sra_run}.log",
    benchmark:
        "logs/benchmarks/SRAdownload/prefetch/{sra_run}.tsv"
    threads: 1
    resources:
        mem=1,
        time=int(config["runtime"]["simplejob"]),
        internet_connection=1,
    conda:
        "%s/sra.yaml" % CONDAENV
    shell:
        " mkdir -p {params.outdir} 2> {log} "
        " ; "
        " prefetch "
        " --output-directory {params.outdir} "
        " -X 999999999 "
        " --progress "
        " --log-level info "
        " {wildcards.sra_run} &>> {log} "
        " ; "
        " vdb-validate {params.outdir}/{wildcards.sra_run}/{wildcards.sra_run}.sra &>> {log} "


rule extract_run:
    input:
        flag=rules.prefetch.output,
    output:
        temp(
            expand(
                SRA_SUBDIR_RUN + "/{{sra_run}}/{{sra_run}}{fraction}.fastq.gz",
                fraction=SRA_read_fractions,
            )
        ),
    params:
        outdir=os.path.abspath(SRA_SUBDIR_RUN + "/{sra_run}"),
        sra_file=SRA_SUBDIR_RUN + "/{sra_run}/{sra_run}.sra",
    log:
        "logs/SRAdownload/extract/{sra_run}.log",
    benchmark:
        "logs/benchmarks/SRAdownload/fasterqdump/{sra_run}.tsv"
    threads: config["simplejob_threads"]
    resources:
        time=int(config["runtime"]["simplejob"]),
        mem=1,  #default 100Mb
    conda:
        "%s/sra.yaml" % CONDAENV
    shell:
        " vdb-validate {params.sra_file} &>> {log} "
        " ; "
        " parallel-fastq-dump "
        " --threads {threads} "
        " --gzip --split-files "
        " --outdir {params.outdir} "
        " --tmpdir {resources.tmpdir} "
        " --skip-technical --split-3 "
        " -s {params.sra_file} &>> {log} "
        " ; "
        " rm -f {params.sra_file} 2>> {log} "


RunTable = None


def get_runids_for_biosample(wildcards):
    global RunTable
    if RunTable is None:

        from atlas.init.parse_sra import load_and_validate_runinfo_table

        RunTable = load_and_validate_runinfo_table("RunInfo.tsv")

    run_ids = RunTable.query(f"BioSample == '{wildcards.sample}'").index.tolist()

    return run_ids


def get_runs_for_biosample(wildcards):

    run_ids = get_runids_for_biosample(wildcards)

    ReadFiles = {}
    for fraction in SRA_read_fractions:

        if fraction == "":
            key = "se"
        else:
            key = fraction

        ReadFiles[key]= expand(
            SRA_SUBDIR_RUN + "/{sra_run}/{sra_run}{fraction}.fastq.gz",
            fraction=fraction,
            sra_run=run_ids,
        )

        return ReadFiles


rule merge_runs_to_sample:
    input:
        unpack(get_runs_for_biosample),
    output:
        expand(
            "SRA/Samples/{{sample}}/{{sample}}{fraction}.fastq.gz",
            fraction=SRA_read_fractions,
        ),
    threads: 1
    run:
        from utils import io

        for i, fraction in enumerate(SRA_read_fractions):

            if fraction == "":
                fraction = "se"
            io.cat_files(input[fraction], output[i])



rule download_sra:
    input:
        expand(
            "SRA/Samples/{sample}/{sample}{fraction}.fastq.gz",
            fraction=SRA_read_fractions,
            sample=SAMPLES,
        ),
