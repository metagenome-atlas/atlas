import os
import re
import sys
from glob import glob
from snakemake.utils import report
import warnings


localrules:
    build_qc_report,
    combine_read_length_stats,
    combine_read_counts,


def get_ribosomal_rna_input(wildcards):
    data_type = config["data_type"]

    clean_reads = expand(
        "{sample}/sequence_quality_control/{sample}_{step}_{fraction}.fastq.gz",
        step="clean",
        fraction=MULTIFILE_FRACTIONS,
        sample=wildcards.sample,
    )
    rrna_reads = expand(
        "{sample}/sequence_quality_control/contaminants/rRNA_{fraction}.fastq.gz",
        fraction=MULTIFILE_FRACTIONS,
        sample=wildcards.sample,
    )

    if data_type == "metagenome" and "rrna" in [
        c.lower() for c in config["contaminant_references"].keys()
    ]:
        return {"clean_reads": clean_reads, "rrna_reads": rrna_reads}
    else:
        return {"clean_reads": clean_reads}


if SKIP_QC:
    PROCESSED_STEPS = ["QC"]

    get_input_fastq = get_quality_controlled_reads

else:
    PROCESSED_STEPS = ["raw"]

    def get_input_fastq(wildcards):
        if not config.get("interleaved_fastqs", False):
            Raw_Headers = ["Reads_raw_" + f for f in RAW_INPUT_FRACTIONS]
        else:
            Raw_Headers = ["Reads_raw_se"]

        # get file
        if not (wildcards.sample in sampleTable.index):
            return expand(
                "Impossible/file/{sample}_{fraction}.fastq.gz",
                sample=wildcards.sample,
                fraction=RAW_INPUT_FRACTIONS,
            )
            logger.debug(
                f"Searched for qc reads for inexisitng sample. wildcards: {wildcards}"
            )
        else:
            return get_files_from_sampleTable(wildcards.sample, Raw_Headers)


rule initialize_qc:
    input:
        unpack(get_input_fastq),
    output:
        temp(
            expand(
                "{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                fraction=RAW_INPUT_FRACTIONS,
                step=PROCESSED_STEPS[-1],
            )
        ),
    priority: 80
    params:
        inputs=lambda wc, input: io_params_for_tadpole(input, "in"),
        interleaved=lambda wc: "t" if config.get("interleaved_fastqs", False) else "f",
        outputs=lambda wc, output: io_params_for_tadpole(output, "out"),
        verifypaired="t" if PAIRED_END else "f",
        extra=config["importqc_params"],
    log:
        "{sample}/logs/QC/init.log",
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads: config.get("simplejob_threads", 1)
    resources:
        mem=config["simplejob_mem"],
        java_mem=int(config["simplejob_mem"] * JAVA_MEM_FRACTION),
    shell:
        "reformat.sh "
        " {params.inputs} "
        " interleaved={params.interleaved} "
        " {params.outputs} "
        " {params.extra} "
        " overwrite=true "
        " verifypaired={params.verifypaired} "
        " threads={threads} "
        " -Xmx{resources.java_mem}G "
        " 2> {log}"


rule get_read_stats:
    input:
        expand(
            "{{sample}}/sequence_quality_control/{{sample}}_{{step}}_{fraction}.fastq.gz",
            fraction=RAW_INPUT_FRACTIONS,
        ),
    output:
        "{sample}/sequence_quality_control/read_stats/{step}.zip",
        read_counts=temp(
            "{sample}/sequence_quality_control/read_stats/{step}_read_counts.tsv"
        ),
    priority: 30
    log:
        "{sample}/logs/QC/read_stats/{step}.log",
    # conda:
    #     "%s/required_packages.yaml" % CONDAENV
    threads: config.get("simplejob_threads", 1)
    resources:
        mem=config["simplejob_mem"],
        java_mem=int(config["simplejob_mem"] * JAVA_MEM_FRACTION),
    params:
        folder=lambda wc, output: os.path.splitext(output[0])[0],
        single_end_file=(
            "{sample}/sequence_quality_control/{sample}_{step}_se.fastq.gz"
        ),
    script:
        "../scripts/get_read_stats.py"


if not SKIP_QC:
    if config.get("deduplicate", True):
        PROCESSED_STEPS.append("deduplicated")

        rule deduplicate_reads:
            input:
                expand(
                    "{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                    step=PROCESSED_STEPS[-2],
                    fraction=RAW_INPUT_FRACTIONS,
                ),
            output:
                temp(
                    expand(
                        "{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                        fraction=RAW_INPUT_FRACTIONS,
                        step=PROCESSED_STEPS[-1],
                    )
                ),
            benchmark:
                "logs/benchmarks/QC/deduplicate/{sample}.txt"
            params:
                inputs=lambda wc, input: io_params_for_tadpole(input, "in"),
                outputs=lambda wc, output: io_params_for_tadpole(output, "out"),
                dupesubs=config.get("duplicates_allow_substitutions"),
                only_optical=("t" if config.get("duplicates_only_optical") else "f"),
            log:
                sterr="{sample}/logs/QC/deduplicate.err",
                stout="{sample}/logs/QC/deduplicate.log",
            conda:
                "%s/required_packages.yaml" % CONDAENV
            threads: config.get("threads", 1)
            resources:
                mem=config["mem"],
                java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
            shell:
                "clumpify.sh "
                " {params.inputs} "
                " {params.outputs} "
                " overwrite=true"
                " dedupe=t "
                " dupesubs={params.dupesubs} "
                " optical={params.only_optical}"
                " threads={threads} "
                " pigz=t unpigz=t "
                " -Xmx{resources.java_mem}G "
                " 2> {log.sterr} "
                " 1> {log.stout} "

    PROCESSED_STEPS.append("filtered")

    rule apply_quality_filter:
        input:
            reads=expand(
                "{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                fraction=RAW_INPUT_FRACTIONS,
                step=PROCESSED_STEPS[-2],
            ),
            adapters=ancient(config["preprocess_adapters"]),
        output:
            reads=temp(
                expand(
                    "{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                    fraction=MULTIFILE_FRACTIONS,
                    step=PROCESSED_STEPS[-1],
                )
            ),
            stats="{sample}/logs/{sample}_quality_filtering_stats.txt",
        benchmark:
            "logs/benchmarks/QC/quality_filter/{sample}.txt"
        params:
            ref=(
                "ref=%s" % config.get("preprocess_adapters")
                if (config["preprocess_adapters"] is not None)
                else ""
            ),
            mink="mink=%d" % config.get("preprocess_adapter_min_k"),
            ktrim="ktrim=%s" % config.get("preprocess_kmer_trim"),
            trimq=config.get("preprocess_minimum_base_quality"),
            hdist="hdist=%d" % config.get("preprocess_allowable_kmer_mismatches"),
            k="k=%d" % config.get("preprocess_reference_kmer_match_length"),
            qtrim=config.get("qtrim"),
            error_correction_pe=(
                "t"
                if PAIRED_END and config["error_correction_overlapping_pairs"]
                else "f"
            ),
            minlength=config["preprocess_minimum_passing_read_length"],
            minbasefrequency=config["preprocess_minimum_base_frequency"],
            # we require the user to reformat to R1 and R2, non-interleaved files
            interleaved="f",
            maxns=config.get("preprocess_max_ns"),
            prealloc=config.get("preallocate_ram"),
            inputs=lambda wc, input: io_params_for_tadpole(input.reads),
            outputs=lambda wc, output: io_params_for_tadpole(
                output.reads, key="out", allow_singletons=False
            ),
        log:
            sterr="{sample}/logs/QC/quality_filter.err",
            stout="{sample}/logs/QC/quality_filter.log",
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads: config.get("threads", 1)
        resources:
            mem=config["mem"],
            java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
        shell:
            " bbduk.sh {params.inputs} "
            " {params.ref} "
            " interleaved={params.interleaved} "
            " {params.outputs} "
            " stats={output.stats} "
            " overwrite=true "
            " qout=33 "
            " trd=t "
            " {params.hdist} "
            " {params.k} "
            " {params.ktrim} "
            " {params.mink} "
            " trimq={params.trimq} "
            " qtrim={params.qtrim} "
            " threads={threads} "
            " minlength={params.minlength} "
            " maxns={params.maxns} "
            " minbasefrequency={params.minbasefrequency} "
            " ecco={params.error_correction_pe} "
            " prealloc={params.prealloc} "
            " pigz=t unpigz=t "
            " -Xmx{resources.java_mem}G "
            " machineout=t "
            " 2> {log.sterr} "
            " 1> {log.stout} "

    # if there are no references, decontamination will be skipped
    if len(config.get("contaminant_references", {}).keys()) > 0:
        PROCESSED_STEPS.append("clean")

        rule build_decontamination_db:
            input:
                ancient(config["contaminant_references"].values()),
            output:
                "ref/genome/1/summary.txt",
            threads: config.get("threads", 1)
            resources:
                mem=config["mem"],
                java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
            log:
                "logs/QC/build_decontamination_db.log",
            conda:
                "%s/required_packages.yaml" % CONDAENV
            params:
                k=config.get("contaminant_kmer_length"),
                refs_in=" ".join(
                    [
                        "ref_%s=%s" % (n, fa)
                        for n, fa in config["contaminant_references"].items()
                    ]
                ),
            shell:
                "bbsplit.sh"
                " -Xmx{resources.java_mem}G "
                " {params.refs_in} "
                " threads={threads}"
                " k={params.k}"
                " local=t "
                " &> {log}"

        rule run_decontamination:
            input:
                reads=expand(
                    "{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                    step=PROCESSED_STEPS[-2],
                    fraction=MULTIFILE_FRACTIONS,
                ),
                db="ref/genome/1/summary.txt",
            output:
                reads=temp(
                    expand(
                        "{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                        fraction=MULTIFILE_FRACTIONS,
                        step=PROCESSED_STEPS[-1],
                    )
                ),
                contaminants=expand(
                    "{{sample}}/sequence_quality_control/contaminants/{db}_{fraction}.fastq.gz",
                    db=list(config["contaminant_references"].keys()),
                    fraction=MULTIFILE_FRACTIONS,
                ),
                stats="{sample}/sequence_quality_control/{sample}_decontamination_reference_stats.txt",
            benchmark:
                "logs/benchmarks/QC/decontamination/{sample}.txt"
            params:
                contaminant_folder=lambda wc, output: os.path.dirname(
                    output.contaminants[0]
                ),
                maxindel=config.get("contaminant_max_indel"),
                minratio=config.get("contaminant_min_ratio"),
                minhits=config.get("contaminant_minimum_hits"),
                ambiguous=config.get("contaminant_ambiguous"),
                k=config.get("contaminant_kmer_length"),
                paired="true" if PAIRED_END else "false",
                inputs=lambda wc, input: io_params_for_tadpole(
                    input.reads, key="in", allow_singletons=False
                ),
                outputs=lambda wc, output: io_params_for_tadpole(
                    output.reads, key="outu", allow_singletons=False
                ),
            log:
                sterr="{sample}/logs/QC/decontamination.err",
                stout="{sample}/logs/QC/decontamination.log",
            conda:
                "../envs/required_packages.yaml"
            threads: config.get("threads", 1)
            resources:
                mem=config["mem"],
                java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
            shell:
                " bbsplit.sh "
                " {params.inputs} "
                " {params.outputs} "
                " basename={params.contaminant_folder}/%_R#.fastq.gz "
                " maxindel={params.maxindel} "
                " minratio={params.minratio} "
                " minhits={params.minhits} "
                " ambiguous={params.ambiguous} "
                " refstats={output.stats} "
                " threads={threads} "
                " k={params.k} "
                " local=t "
                " machineout=t "
                " pigz=t unpigz=t ziplevel=9 "
                " -Xmx{resources.java_mem}G "
                " 1> {log.stout} "
                " 2> {log.sterr} "

    PROCESSED_STEPS.append("QC")

    localrules:
        qcreads,

    rule qcreads:
        input:
            unpack(get_ribosomal_rna_input),
        output:
            temp(
                expand(
                    "{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                    fraction=MULTIFILE_FRACTIONS,
                    step=PROCESSED_STEPS[-1],
                )
            ),
        threads: 1
        run:
            import shutil
            import pandas as pd

            for i in range(len(MULTIFILE_FRACTIONS)):
                with open(output[i], "wb") as outFile:
                    with open(input.clean_reads[i], "rb") as infile1:
                        shutil.copyfileobj(infile1, outFile)
                        if hasattr(input, "rrna_reads"):
                            with open(input.rrna_reads[i], "rb") as infile2:
                                shutil.copyfileobj(infile2, outFile)



rule copy_qc_reads:
    input:
        reads=expand(
            "{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
            fraction=MULTIFILE_FRACTIONS,
            step="QC",
        ),
    output:
        reads=expand(
            "QC/reads/{{sample}}_{fraction}.fastq.gz",
            fraction=MULTIFILE_FRACTIONS,
        ),
    run:
        for i, f in enumerate(input.reads):
            shutil.copy(f, output.reads[i])


#### STATS


if PAIRED_END:

    rule calculate_insert_size:
        input:
            get_quality_controlled_reads,
        output:
            ihist=(
                "{sample}/sequence_quality_control/read_stats/QC_insert_size_hist.txt"
            ),
            read_length=(
                "{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt"
            ),
        threads: config.get("simplejob_threads", 1)
        resources:
            mem=config["mem"],
            java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
        conda:
            "../envs/required_packages.yaml"
        log:
            "{sample}/logs/QC/stats/calculate_insert_size.log",
        params:
            kmer=config.get("merging_k"),
            extend2=config.get("merging_extend2"),
            flags="loose ecct",
            minprob=config.get("bbmerge_minprob", "0.8"),
            inputs=lambda wc, input: io_params_for_tadpole(input),
        shell:
            " bbmerge.sh "
            " -Xmx{resources.java_mem}G "
            " threads={threads} "
            " {params.inputs} "
            " {params.flags} k={params.kmer} "
            " extend2={params.extend2} "
            " ihist={output.ihist} merge=f "
            " mininsert0=35 minoverlap0=8 "
            " prealloc=t prefilter=t "
            " minprob={params.minprob} 2> {log} \n  "
            """
            readlength.sh {params.inputs} out={output.read_length} 2>> {log}
            """

else:

    rule calculate_read_length_hist:
        input:
            get_quality_controlled_reads,
        output:
            read_length=(
                "{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt"
            ),
        params:
            kmer=config.get("merging_k"),
        threads: config.get("simplejob_threads")
        resources:
            mem=config["simplejob_mem"],
        conda:
            "../envs/required_packages.yaml"
        log:
            "{sample}/logs/QC/stats/calculate_read_length.log",
        shell:
            """
            readlength.sh in={input[0]} out={output.read_length} 2> {log}
            """


rule combine_read_length_stats:
    input:
        expand(
            "{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt",
            sample=SAMPLES,
        ),
    output:
        "stats/read_length_stats.tsv",
    run:
        import pandas as pd
        import os
        from utils.parsers_bbmap import parse_comments

        stats = pd.DataFrame()

        for length_file in input:
            sample = length_file.split(os.path.sep)[0]
            data = parse_comments(length_file)
            data = pd.Series(data)[
                ["Reads", "Bases", "Max", "Min", "Avg", "Median", "Mode", "Std_Dev"]
            ]
            stats[sample] = data

        stats.to_csv(output[0], sep="\t")



# rule combine_cardinality:
#     input:
#         expand("{sample}/sequence_quality_control/read_stats/QC_cardinality.txt",sample=SAMPLES),
#     output:
#         'stats/cardinality.tsv'
#     run:
#         import pandas as pd
#         import os

#         stats= pd.Series()

#         for file in input:
#             sample= file.split(os.path.sep)[0]
#             with open(file) as f:
#                 cardinality= int(f.read().strip())

#             stats.loc[sample]=cardinality

#         stats.to_csv(output[0],sep='"t')


if PAIRED_END:

    localrules:
        combine_insert_stats,

    rule combine_insert_stats:
        input:
            expand(
                "{sample}/sequence_quality_control/read_stats/QC_insert_size_hist.txt",
                sample=SAMPLES,
            ),
        output:
            "stats/insert_stats.tsv",
        run:
            import pandas as pd
            import os
            from utils.parsers_bbmap import parse_comments

            stats = pd.DataFrame()

            for insert_file in input:
                sample = insert_file.split(os.path.sep)[0]
                data = parse_comments(insert_file)
                data = pd.Series(data)[
                    ["Mean", "Median", "Mode", "STDev", "PercentOfPairs"]
                ]
                stats[sample] = data

            stats.T.to_csv(output[0], sep="\t")



localrules:
    combine_read_counts,
    write_read_counts,


rule write_read_counts:
    input:
        read_count_files=expand(
            "{{sample}}/sequence_quality_control/read_stats/{step}_read_counts.tsv",
            step=PROCESSED_STEPS,
        ),
    output:
        read_stats="{sample}/sequence_quality_control/read_stats/read_counts.tsv",
    run:
        from utils.io import pandas_concat

        pandas_concat(
            list(input.read_count_files),
            output.read_stats,
            sep="\t",
            index_col=[0, 1],
            axis=0,
        )


rule combine_read_counts:
    input:
        expand(
            "{sample}/sequence_quality_control/read_stats/read_counts.tsv",
            sample=SAMPLES,
        ),
    output:
        "stats/read_counts.tsv",
    run:
        from utils.io import pandas_concat

        pandas_concat(list(input), output[0], sep="\t", index_col=[0, 1], axis=0)


rule finalize_sample_qc:
    input:
        reads=expand(
            "QC/reads/{{sample}}_{fraction}.fastq.gz",
            fraction=MULTIFILE_FRACTIONS,
        ),
        #quality_filtering_stats = "{sample}/logs/{sample}_quality_filtering_stats.txt",
        reads_stats_zip=expand(
            "{{sample}}/sequence_quality_control/read_stats/{step}.zip",
            step=PROCESSED_STEPS,
        ),
        read_length_hist=(
            "{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt"
        ),
    output:
        flag=touch("{sample}/sequence_quality_control/finished_QC"),


rule build_qc_report:
    input:
        zipfiles_QC=expand(
            "{sample}/sequence_quality_control/read_stats/QC.zip", sample=SAMPLES
        ),
        read_counts="stats/read_counts.tsv",
        read_length_stats=(
            ["stats/read_length_stats.tsv", "stats/insert_stats.tsv"]
            if PAIRED_END
            else "stats/read_length_stats.tsv"
        ),
    output:
        report="reports/QC_report.html",
    log:
        "logs/QC/report.log",
    params:
        min_quality=config["preprocess_minimum_base_quality"],
        samples=SAMPLES,
    conda:
        "../envs/report.yaml"
    script:
        "../report/qc_report.py"
