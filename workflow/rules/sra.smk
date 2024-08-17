wildcard_constraints:
    sra_run="[S,E,D]RR[0-9]+",


localrules:
    kingfisher_get,
    merge_runs_to_sample


SRA_read_fractions = ["_1", "_2"] if PAIRED_END else [""]
SRA_SUBDIR_RUN = Path("SRA/Runs")


RunTable = None
def load_runtable():
    global RunTable
    if RunTable is None:
        from atlas.init import parse_sra
        RunTable = parse_sra.load_and_validate_runinfo_table()
    return RunTable

def get_run_ids_for_sample(wildcards):
    RunTable = load_runtable()
    from atlas.init import parse_sra
    return parse_sra.get_run_ids_for_sample(RunTable, wildcards.sample)


rule kingfisher_get:
    output:
        #dir = temp(directory("Reads/tmp/runs/{sample}")),
        flag = temp(touch("Reads/tmp/flags/{sample}.downloaded")),
    params:
        run_ids = get_run_ids_for_sample,
        download_methods="ena-ascp ena-ftp prefetch",
        output_dir= lambda wc: SRA_SUBDIR_RUN / wc.sample,
    log:
        Path("log/download_reads/download/{sample}.log").resolve(),
    threads: config['threads'],
    resources:
        mem_mb=config['mem']*1000,
        time_min=config["runtime"]["long"]*60,
        ncbi_connection=1
    conda:
        "../envs/kingfisher.yaml"
    shell:
        " mkdir -p {params.output_dir} ; "
        " cd {params.output_dir} "
        " ; "
        "kingfisher get --run-identifiers {params.run_ids} "
        " --download-threads 2 --extraction-threads {threads} "
        " --hide-download-progress "
        " --output-format-possibilities 'fastq.gz' "
        " --force --check-md5sums "
        " --download-methods {params.download_methods} "
        " -f fastq.gz &> {log}"



def get_run_fastq_for_sample(run_ids):

    ReadFiles = {}
    for fraction in SRA_read_fractions:
        if fraction == "":
            key = "se"
        else:
            key = fraction

        ReadFiles[key] = expand(
            str(SRA_SUBDIR_RUN / "{sample}/{sra_run}{fraction}.fastq.gz"),
            fraction=fraction,
            sra_run=run_ids,
        )

    return ReadFiles


rule merge_runs_to_sample:
    input:
        flag= "Reads/tmp/flags/{sample}.downloaded"
    output:
        expand(
            "SRA/Samples/{{sample}}/{{sample}}{fraction}.fastq.gz",
            fraction=SRA_read_fractions,
        ),
    params:
        run_ids = get_run_ids_for_sample,
    threads: 1
    run:

        # print(list( (SRA_SUBDIR_RUN / wildcards.sample).iterdir()))

        from utils import io


        for i, fraction in enumerate(SRA_read_fractions):

            run_fastqs= expand(
                SRA_SUBDIR_RUN / "{sample}/{sra_run}{fraction}.fastq.gz",
                fraction=fraction,
                sra_run=params.run_ids,
                sample=wildcards.sample
            )

            assert all([Path(f).exists() for f in run_fastqs])," Not all fastq files exist. Expected: %s" % run_fastqs

            io.cat_files(run_fastqs, output[i])







rule download_sra:
    input:
        expand(
            "SRA/Samples/{sample}/{sample}{fraction}.fastq.gz",
            fraction=SRA_read_fractions,
            sample=SAMPLES,
        ),
