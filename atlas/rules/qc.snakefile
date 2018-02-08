import os
import re
import sys
from glob import glob
from snakemake.utils import report
import warnings


localrules: postprocess_after_decontamination, initialize_checkm, \
            finalize_QC, QC_report, combine_read_length_stats, \
            combine_insert_stats, combine_read_counts


def get_line_from_file(file, line_start):
    "function to extract line after string name, eg in datafile with average in comments #Avg: 123"

    with open(file) as f:
        for line in f:
            if line.startswith(line_start):
                break

    if not line.startswith(line_start):
        raise Exception("Didn't find {name} in file ({file}):\n\n".format(name=name,file=file))
    else:
        return line[len(line_start):].strip()


def parse_comments(file, comment='#',sep='\t',expect_one_value=True):
    "parse comments at begin of file #Avg: 123"
    Parsed = {}
    with open(file) as f:
        line = f.readline()

        while line.startswith(comment):
            line_values = line[1:].strip().split(sep)
            name = line_values[0]
            if name[-1] == ':':
                name = name[:-1]
            values = line_values[1:]

            if len(values) == 1:
                Parsed[name]=values[0]
            elif not expect_one_value:
                Parsed[name]=values

            line= f.readline()

    if len(Parsed)==0:
        raise Exception("Couldn't parse values from file {} with comment char {} and sep '{}' ".format(file,comment,sep))
    return Parsed


def get_ribosomal_rna_input(wildcards):
    data_type = config["samples"][wildcards.sample].get("type", "metagenome").lower()

    clean_reads = expand("{sample}/sequence_quality_control/{sample}_{step}_{fraction}.fastq.gz",step='clean', fraction=MULTIFILE_FRACTIONS,sample=wildcards.sample)
    rrna_reads = expand("{sample}/sequence_quality_control/contaminants/rRNA_{fraction}.fastq.gz",fraction=MULTIFILE_FRACTIONS,sample=wildcards.sample)

    if data_type == "metagenome" and 'rrna' in [c.lower() for c in config["contaminant_references"].keys()]:
        return {'clean_reads':clean_reads, 'rrna_reads':rrna_reads}
    else:
        return {'clean_reads':clean_reads}


def get_finalize_qc_input(wildcards):
    inputs = get_quality_controlled_reads(wildcards)
    try:
        # FIXME: 'decontaminated_seqs': ['{sample}/sequence_quality_control/contaminants/PhiX_R1.fastq.gz'...
        # inputs["decontaminated_seqs"] = rules.decontamination.output.contaminants
        inputs["decontaminated_stats"] = "{sample}/sequence_quality_control/{sample}_decontamination_reference_stats.txt".format(sample=wildcards.sample)
    except AttributeError:
        pass
    return inputs


PROCESSED_STEPS = ['raw']


# controls files and deinterlevves them, for the pipeline all files have the same format
rule init_QC:
    input:
        lambda wc: config["samples"][wc.sample]["fastq"]
    output:
        temp(expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                    fraction=RAW_INPUT_FRACTIONS, step=PROCESSED_STEPS[-1]))
    priority: 80
    params:
        inputs = lambda wc: "in=%s" % config["samples"][wc.sample]["fastq"][0] if len(config["samples"][wc.sample]["fastq"]) == 1 else "in=%s in2=%s" % tuple(config["samples"][wc.sample]["fastq"]),
        interleaved = lambda wc: "t" if config["samples"][wc.sample].get("paired", True) and len(config["samples"][wc.sample]["fastq"]) == 1 else "f",
        outputs = lambda wc,output: "out1={0} out2={1}".format(*output) if PAIRED_END else "out={0}".format(*output),
        verifypaired = "t" if PAIRED_END else "f"
    log:
        "{sample}/logs/{sample}_init.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", JAVA_MEM),
        java_mem = int(int(config.get("java_mem", JAVA_MEM) * JAVA_MEM_FRACTION))
    shell:
        """reformat.sh {params.inputs} \
        interleaved={params.interleaved} \
        {params.outputs} \
        qout=33 \
        overwrite=true \
        verifypaired={params.verifypaired} \
        addslash=t \
        trimreaddescription=t \
        threads={threads} \
        -Xmx{resources.java_mem}G 2> {log}
        """


rule read_stats:
    # TODO: remove run block in favor of script or alternate cli
    # see http://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#snakefiles-external-scripts
    input:
        expand("{{sample}}/sequence_quality_control/{{sample}}_{{step}}_{fraction}.fastq.gz",
            fraction=RAW_INPUT_FRACTIONS)
    output:
        "{sample}/sequence_quality_control/read_stats/{step}.zip",
        read_counts = temp("{sample}/sequence_quality_control/read_stats/{step}_read_counts.tsv"),
    priority:
        30
    log:
        "{sample}/logs/read_stats.log"
    # conda:
    #     "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", JAVA_MEM),
        java_mem = int(config.get("java_mem", JAVA_MEM) * JAVA_MEM_FRACTION)
    params:
        folder = lambda wc, output: os.path.splitext(output[0])[0],
        single_end_file = "{sample}/sequence_quality_control/{sample}_{step}_se.fastq.gz"
    run:
        import datetime
        import shutil

        timestamp = datetime.datetime.now().strftime('%Y-%m-%d-%X')

        def get_read_stats(fraction, params_in):

            subfolder = os.path.join(params.folder, fraction)
            tmp_file=os.path.join(subfolder,"read_stats.tmp")
            shell("""
                    mkdir -p {subfolder}

                    reformat.sh {params_in} \
                    bhist={subfolder}/base_hist.txt \
                    qhist={subfolder}/quality_by_pos.txt \
                    lhist={subfolder}/readlength.txt \
                    gchist={subfolder}/gc_hist.txt \
                    gcbins=auto \
                    bqhist={subfolder}/boxplot_quality.txt \
                    threads={threads} \
                    overwrite=true \
                    -Xmx{mem}G \
                    2> >(tee -a {log} {tmp_file} )
                 """.format(subfolder=subfolder, params_in=params_in, log=log,
                            threads=threads, mem=resources.java_mem,tmp_file=tmp_file))
            content = open(tmp_file).read()
            pos = content.find('Input:')
            if pos == -1:
                raise Exception("Didn't find read number in file:\n\n" + content)
            else:

                content[pos:].split()[1:4]
                        # Input:    123 reads   1234 bases
                n_reads, _, n_bases = content[pos:].split()[1:4]

                os.remove(tmp_file)
            return int(n_reads), int(n_bases)


        if PAIRED_END:
            n_reads_pe, n_bases_pe = get_read_stats('pe', "in1={0} in2={1}".format(*input))
            n_reads_pe = n_reads_pe / 2
            headers = ['Sample', 'Step', 'Total_Reads', 'Total_Bases',
                       'Reads_pe', 'Bases_pe', 'Reads_se', 'Bases_se',
                       'Timestamp']

            if os.path.exists(params.single_end_file):
                n_reads_se, n_bases_se = get_read_stats('se', "in=" + params.single_end_file)
            else:
                n_reads_se, n_bases_se = 0, 0

            values=[n_reads_pe + n_reads_se, n_bases_pe + n_bases_se,
                    n_reads_pe, n_bases_pe,
                    n_reads_se, n_bases_se]
        else:
            headers= ['Sample', 'Step', 'Total_Reads', 'Total_Bases', 'Reads',
                      'Bases', 'Timestamp']
            values = 2 * get_read_stats('', "in=" + input[0])

        with open(output.read_counts, 'w') as f:
            f.write('\t'.join(headers) + '\n')
            f.write('\t'.join([wildcards.sample, wildcards.step] + [str(v) for v in values] + [timestamp]) + '\n')

        shutil.make_archive(params.folder, 'zip', params.folder)
        shutil.rmtree(params.folder)


if config.get('deduplicate', True):
    PROCESSED_STEPS.append("deduplicated")
    rule deduplicate:
        input:
            expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                step=PROCESSED_STEPS[-2], fraction=RAW_INPUT_FRACTIONS)
        output:
            temp(expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                fraction=RAW_INPUT_FRACTIONS, step=PROCESSED_STEPS[-1]))
        benchmark:
            "logs/benchmarks/deduplicate/{sample}.txt"
        params:
            inputs = lambda wc, input: "in=%s in2=%s" % (input[0], input[1]) if PAIRED_END else "in=%s" % input[0],
            outputs = lambda wc,output: "out1={0} out2={1}".format(*output) if PAIRED_END else "out={0}".format(*output),
            dupesubs = config.get('duplicates_allow_substitutions', DUPLICATES_ALLOW_SUBSTITUTIONS),
            only_optical = 't' if config.get('duplicates_only_optical', DUPLICATES_ONLY_OPTICAL) else 'f'
        log:
            "{sample}/logs/{sample}_deduplicate.log"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        resources:
            mem = config.get("java_mem", JAVA_MEM),
            java_mem = int(config.get("java_mem", JAVA_MEM) * JAVA_MEM_FRACTION)
        shell:
            """
            clumpify.sh \
                {params.inputs} \
                {params.outputs} \
                overwrite=true\
                dedupe=t \
                dupesubs={params.dupesubs} \
                optical={params.only_optical}\
                threads={threads} \
                -Xmx{resources.java_mem}G 2> {log}
            """


PROCESSED_STEPS.append("filtered")


rule quality_filter:
    input:
        expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
            fraction=RAW_INPUT_FRACTIONS, step=PROCESSED_STEPS[-2])
    output:
        temp(expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
            fraction=MULTIFILE_FRACTIONS, step=PROCESSED_STEPS[-1])),
        stats = "{sample}/logs/{sample}_quality_filtering_stats.txt"
    benchmark:
        "logs/benchmarks/quality_filter/{sample}.txt"
    params:
        ref = "ref=%s" % config.get("preprocess_adapters") if config.get("preprocess_adapters") else "",
        mink = "" if not config.get("preprocess_adapters") else "mink=%d" % config.get("preprocess_adapter_min_k", PREPROCESS_ADAPTER_MIN_K),
        ktrim = "" if not config.get("preprocess_adapters") else "ktrim=%s" % config.get("preprocess_kmer_trim", PREPROCESS_KMER_TRIM),
        trimq = config.get("preprocess_minimum_base_quality", PREPROCESS_MINIMUM_BASE_QUALITY),
        hdist = "" if not config.get("preprocess_adapters") else "hdist=%d" % config.get("preprocess_allowable_kmer_mismatches", PREPROCESS_ALLOWABLE_KMER_MISMATCHES),
        k = "" if not config.get("preprocess_adapters") else "k=%d" % config.get("preprocess_reference_kmer_match_length", PREPROCESS_REFERENCE_KMER_MATCH_LENGTH),
        qtrim = config.get("qtrim", QTRIM),
        error_correction_pe = "t" if PAIRED_END and config.get('error_correction_overlapping_pairs',True) else "f",
        minlength = config.get("preprocess_minimum_passing_read_length", PREPROCESS_MINIMUM_PASSING_READ_LENGTH),
        minbasefrequency = config.get("preprocess_minimum_base_frequency", PREPROCESS_MINIMUM_BASE_FREQUENCY),
        # we require the user to reformat to R1 and R2, non-interleaved files
        interleaved = "f",
        maxns = config.get("preprocess_max_ns", PREPROCESS_MAX_NS),
        prealloc = config.get("preallocate_ram", PREALLOCATE_RAM),
        inputs = lambda wc, input:"in={0} in2={1}".format(*input) if PAIRED_END else "in={0}".format(*input),
        outputs = lambda wc, output:"out={0} out2={1} outs={2}".format(*output) if PAIRED_END else "out={0}".format(*output)
    log:
        "{sample}/logs/{sample}_quality_filter.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", JAVA_MEM),
        java_mem = int(config.get("java_mem", JAVA_MEM) * JAVA_MEM_FRACTION)
    shell:
        """
        bbduk.sh {params.inputs} \
            {params.ref} \
            interleaved={params.interleaved} \
            {params.outputs} \
            stats={output.stats} \
            overwrite=true \
            qout=33 \
            trd=t \
            {params.hdist} \
            {params.k} \
            {params.ktrim} \
            {params.mink} \
            trimq={params.trimq} \
            qtrim={params.qtrim} \
            threads={threads} \
            minlength={params.minlength} \
            maxns={params.maxns} \
            minbasefrequency={params.minbasefrequency} \
            ecco={params.error_correction_pe} \
            prealloc={params.prealloc} \
            -Xmx{resources.java_mem}G 2> {log}
        """


# if there are no references, decontamination will be skipped
if len(config.get("contaminant_references", {}).keys()) > 0:
    PROCESSED_STEPS.append("clean")

    rule build_decontamination_db:
        output:
            "ref/genome/1/summary.txt"
        threads:
            config.get("threads", 1)
        resources:
            mem = config.get("java_mem", JAVA_MEM),
            java_mem = int(config.get("java_mem", JAVA_MEM) * JAVA_MEM_FRACTION)
        log:
            "logs/build_decontamination_db.log"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        params:
            k = config.get("contaminant_kmer_length", CONTAMINANT_KMER_LENGTH),
            refs_in = " ".join(["ref_%s=%s" % (n, fa) for n, fa in config["contaminant_references"].items()]),
        shell:
            """bbsplit.sh -Xmx{resources.java_mem}G {params.refs_in} threads={threads} k={params.k} local=t 2> {log}"""


    rule decontamination:
        input:
            expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                step=PROCESSED_STEPS[-2], fraction=MULTIFILE_FRACTIONS),
            db = "ref/genome/1/summary.txt"
        output:
            temp(expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                fraction=MULTIFILE_FRACTIONS, step=PROCESSED_STEPS[-1])),
            contaminants = expand("{{sample}}/sequence_quality_control/contaminants/{db}_{fraction}.fastq.gz",
                    db=list(config["contaminant_references"].keys()),
                    fraction=MULTIFILE_FRACTIONS),
            stats = "{sample}/sequence_quality_control/{sample}_decontamination_reference_stats.txt"
        benchmark:
            "logs/benchmarks/decontamination/{sample}.txt"
        params:
            contaminant_folder = lambda wc, output: os.path.dirname(output.contaminants[0]),
            maxindel = config.get("contaminant_max_indel", CONTAMINANT_MAX_INDEL),
            minratio = config.get("contaminant_min_ratio", CONTAMINANT_MIN_RATIO),
            minhits = config.get("contaminant_minimum_hits", CONTAMINANT_MINIMUM_HITS),
            ambiguous = config.get("contaminant_ambiguous", CONTAMINANT_AMBIGUOUS),
            k = config.get("contaminant_kmer_length", CONTAMINANT_KMER_LENGTH),
            paired= "true" if PAIRED_END else "false",
            input_single= lambda wc,input: input[2] if PAIRED_END else input[0],
            output_single= lambda wc,output: output[2] if PAIRED_END else output[0]
        log:
            "{sample}/logs/{sample}_decontamination.log"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        resources:
            mem = config.get("java_mem", JAVA_MEM),
            java_mem = int(config.get("java_mem", JAVA_MEM) * JAVA_MEM_FRACTION)
        shell:
            """
            if [ "{params.paired}" = true ] ; then
                bbsplit.sh in1={input[0]} in2={input[1]} \
                    outu1={output[0]} outu2={output[1]} \
                    basename="{params.contaminant_folder}/%_R#.fastq.gz" \
                    maxindel={params.maxindel} minratio={params.minratio} \
                    minhits={params.minhits} ambiguous={params.ambiguous} refstats={output.stats}\
                    threads={threads} k={params.k} local=t \
                    -Xmx{resources.java_mem}G 2> {log}
            fi

            bbsplit.sh in={params.input_single}  \
                outu={params.output_single} \
                basename="{params.contaminant_folder}/%_se.fastq.gz" \
                maxindel={params.maxindel} minratio={params.minratio} \
                minhits={params.minhits} ambiguous={params.ambiguous} refstats={output.stats} append \
                interleaved=f threads={threads} k={params.k} local=t \
                -Xmx{resources.java_mem}G 2>> {log}
            """


PROCESSED_STEPS.append("QC")


rule postprocess_after_decontamination:
    input:
        unpack(get_ribosomal_rna_input)
    output:
        expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",step=PROCESSED_STEPS[-1],fraction=MULTIFILE_FRACTIONS)
    threads:
        1
    run:
        import shutil
        for i in range(len(MULTIFILE_FRACTIONS)):
            with open(output[i], 'wb') as outFile:
                with open(input.clean_reads[i], 'rb') as infile1:
                    shutil.copyfileobj(infile1, outFile)
                    if hasattr(input, 'rrna_reads'):
                        with open(input.rrna_reads[i], 'rb') as infile2:
                            shutil.copyfileobj(infile2, outFile)

if PAIRED_END:
    rule calculate_insert_size:
        input:
            unpack(get_quality_controlled_reads)
        output:
            ihist = "{sample}/sequence_quality_control/read_stats/QC_insert_size_hist.txt",
            read_length= "{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt"
        threads:
            config.get("threads", 1)
        resources:
            mem = config.get("java_mem", JAVA_MEM),
            java_mem = int(config.get("java_mem", JAVA_MEM) * JAVA_MEM_FRACTION)
        conda:
            "%s/required_packages.yaml" % CONDAENV
        log:
            "{sample}/logs/{sample}_calculate_insert_size.log"
        benchmark:
            "logs/benchmarks/merge_pairs/{sample}_insert_size.txt"
        params:
            kmer = config.get("merging_k", MERGING_K),
            extend2 = config.get("merging_extend2", MERGING_EXTEND2),
            flags = 'loose ecct',
            minprob = config.get("bbmerge_minprob", "0.8")
        shell:
            """
            bbmerge.sh -Xmx{resources.java_mem}G threads={threads} \
                in1={input.R1} in2={input.R2} \
                {params.flags} k={params.kmer} \
                extend2={params.extend2} \
                ihist={output.ihist} merge=f \
                mininsert0=35 minoverlap0=8 \
                prealloc=t prefilter=t \
                minprob={params.minprob} 2> >(tee {log})

            readlength.sh in={input.R1} in2={input.R2} out={output.read_length} 2> >(tee {log})
            """


else:
    rule calculate_read_length_hist:
        input:
            unpack(get_quality_controlled_reads)
        output:
            read_length= "{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt",
        params:
            kmer = config.get("merging_k", MERGING_K),
        threads:
            config.get("threads", 1)
        resources:
            mem = config.get("java_mem", JAVA_MEM),
            java_mem = int(config.get("java_mem", JAVA_MEM) * JAVA_MEM_FRACTION)
        conda:
            "%s/required_packages.yaml" % CONDAENV
        log:
            "{sample}/logs/{sample}_calculate_read_length.log"
        shell:
            """
            readlength.sh in={input.se} out={output.read_length} 2> >(tee {log})
            """


rule combine_read_length_stats:
    input:
        expand("{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt", sample=SAMPLES),
    output:
        'stats/read_length_stats.tsv'
    run:
        import pandas as pd
        import os

        stats = pd.DataFrame()

        for length_file in input:
            sample = length_file.split(os.path.sep)[0]
            data = parse_comments(length_file)
            data = pd.Series(data)[['Reads', 'Bases', 'Max', 'Min', 'Avg', 'Median', 'Mode', 'Std_Dev']]
            stats[sample] = data

        stats.to_csv(output[0], sep='\t')


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

#         stats.to_csv(output[0],sep='\t')


if PAIRED_END:
    rule combine_insert_stats:
        input:
            expand("{sample}/sequence_quality_control/read_stats/QC_insert_size_hist.txt",sample=SAMPLES),
        output:
            'stats/insert_stats.tsv'
        run:
            import pandas as pd
            import os

            stats = pd.DataFrame()

            for insert_file in input:
                sample = insert_file.split(os.path.sep)[0]
                data = parse_comments(insert_file)
                data = pd.Series(data)[['Mean', 'Median', 'Mode', 'STDev', 'PercentOfPairs']]
                stats[sample] = data

            stats.T.to_csv(output[0], sep='\t')


rule combine_read_counts:
    input:
        expand("{sample}/sequence_quality_control/read_stats/read_counts.tsv",sample=SAMPLES)
    output:
        "stats/read_counts.tsv"
    run:
        import pandas as pd

        stats = pd.DataFrame()

        for f in input:
            d = pd.read_table(f, index_col=[0, 1])
            stats = stats.append(d)

        stats.to_csv(output[0], sep='\t')


rule finalize_QC:
    input:
        unpack(get_finalize_qc_input),
        quality_filtering_stats = "{sample}/logs/{sample}_quality_filtering_stats.txt",
        reads_stats_zip = expand("{{sample}}/sequence_quality_control/read_stats/{step}.zip", step=PROCESSED_STEPS),
        read_count_files = expand("{{sample}}/sequence_quality_control/read_stats/{step}_read_counts.tsv", step=PROCESSED_STEPS),
        read_length_hist = "{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt"
    output:
        touch("{sample}/sequence_quality_control/finished_QC"),
        read_stats = "{sample}/sequence_quality_control/read_stats/read_counts.tsv"
    run:
        import pandas as pd

        all_read_counts = pd.DataFrame()
        for read_stats_file in input.read_count_files:
            d = pd.read_table(read_stats_file, index_col=[0, 1])
            all_read_counts = all_read_counts.append(d)
        all_read_counts.to_csv(output.read_stats, sep='\t')
        print("Finished QC for sample {sample}\n".format(**wildcards))


rule QC_report:
    input:
        expand("{sample}/sequence_quality_control/finished_QC", sample=SAMPLES),
        read_counts= "stats/read_counts.tsv",
        read_length_stats= ['stats/insert_stats.tsv','stats/read_length_stats.tsv'] if PAIRED_END else 'stats/read_length_stats.tsv'
    output:
        touch("finished_QC"),
        report = "reports/QC_report.html"
    conda:
        "%s/report.yaml" % CONDAENV
    script:
        "../report/qc_report.py"


# aggregate stats reports ...
