import os
import re
import sys
from glob import glob
from snakemake.utils import report
import warnings


localrules: postprocess_after_decontamination,initialize_checkm,finalize_QC,QC_report

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


processed_steps=['raw']

# controls files and deinterlevves them, for the pipeline all files have the same format

rule init_QC:
    input:
        lambda wc: config["samples"][wc.sample]["fastq"]
    output:
        temp(expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
            fraction=raw_input_fractions,step=processed_steps[-1]))
    priority: 80
    params:
        inputs = lambda wc: "in=%s" % config["samples"][wc.sample]["fastq"][0] if len(config["samples"][wc.sample]["fastq"]) == 1 else "in=%s in2=%s" % tuple(config["samples"][wc.sample]["fastq"]),
        interleaved = lambda wc: "t" if config["samples"][wc.sample].get("paired", True) and len(config["samples"][wc.sample]["fastq"]) == 1 else "f",
        outputs = lambda wc,output: "out1={0} out2={1}".format(*output) if paired_end else "out={0}".format(*output),
        verifypaired="t" if paired_end else "f"
    log:
        "{sample}/logs/{sample}_init.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", 5)
    shell:
        """{SHPFXM} reformat.sh {params.inputs} \
        interleaved={params.interleaved} \
        {params.outputs} \
        qout=33 \
        overwrite=true \
        verifypaired={params.verifypaired} \
        addslash=t \
        trimreaddescription=t \
        threads={threads} \
        -Xmx{resources.mem}G 2> {log}
        """

rule read_stats:
    # TODO: remove run block in favor of script or alternate cli
    # see http://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#snakefiles-external-scripts
    input:
        expand("{{sample}}/sequence_quality_control/{{sample}}_{{step}}_{fraction}.fastq.gz",
            fraction=raw_input_fractions)
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
        mem = config.get("java_mem", JAVA_MEM)
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
                            threads=threads, mem=resources.mem,tmp_file=tmp_file))
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


        if paired_end:
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
    processed_steps += ['deduplicated']
    rule deduplicate:
        input:
            expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                step=processed_steps[-2], fraction=raw_input_fractions)
        output:
            temp(expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                fraction=raw_input_fractions, step=processed_steps[-1]))
        benchmark:
            "logs/benchmarks/deduplicate/{sample}.txt"
        params:
            inputs = lambda wc, input: "in=%s in2=%s" % (input[0], input[1]) if paired_end else "in=%s" % input[0],
            outputs = lambda wc,output: "out1={0} out2={1}".format(*output) if paired_end else "out={0}".format(*output),
            dupesubs= config.get('duplicates_allow_substitutions', DUPLICATES_ALLOW_SUBSTITUTIONS),
            only_optical = 't' if config.get('duplicates_only_optical',DUPLICATES_ONLY_OPTICAL) else 'f'
        log:
            "{sample}/logs/{sample}_deduplicate.log"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        resources:
            mem = config.get("java_mem", JAVA_MEM)
        shell:
            """

            {SHPFXM} clumpify.sh \
            {params.inputs} \
            {params.outputs} \
            overwrite=true\
            dedupe=t \
            dupesubs={params.dupesubs} \
            optical={params.only_optical}\
            threads={threads} \
            -Xmx{resources.mem}G 2> {log}

            """



processed_steps += ['filtered']

rule quality_filter:
    input:
        expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
            fraction=raw_input_fractions, step=processed_steps[-2])
    output:
        temp(expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
            fraction=multifile_fractions, step=processed_steps[-1])),
        stats = "{sample}/logs/{sample}_quality_filtering_stats.txt"
    benchmark:
        "logs/benchmarks/quality_filter/{sample}.txt"
    params:
        lref = "lref=%s" % config.get("preprocess_adapters") if config.get("preprocess_adapters") else "",
        rref = "rref=%s" % config.get("preprocess_adapters") if config.get("preprocess_adapters") else "",
        mink = "" if not config.get("preprocess_adapters") else "mink=%d" % config.get("preprocess_adapter_min_k", 8),
        trimq = config.get("preprocess_minimum_base_quality", PREPROCESS_MINIMUM_BASE_QUALITY),
        hdist = "" if not config.get("preprocess_adapters") else "hdist=%d" % config.get("preprocess_allowable_kmer_mismatches", PREPROCESS_ALLOWABLE_KMER_MISMATCHES),
        k = "" if not config.get("preprocess_adapters") else "k=%d" % config.get("preprocess_reference_kmer_match_length", PREPROCESS_REFERENCE_KMER_MATCH_LENGTH),
        qtrim = config.get("qtrim", QTRIM),
        error_correction_pe= "t" if paired_end and config.get('error_correction_overlapping_pairs',True) else "f",
        minlength = config.get("preprocess_minimum_passing_read_length", PREPROCESS_MINIMUM_PASSING_READ_LENGTH),
        minbasefrequency = config.get("preprocess_minimum_base_frequency", PREPROCESS_MINIMUM_BASE_FREQUENCY),
        interleaved = "t" if paired_end else "f",
        inputs= lambda wc,input:"in1={0} in2={1}".format(*input) if paired_end else "in={0}".format(*input),
        outputs=  lambda wc,output:"out1={0} out2={1} outs={2}".format(*output) if paired_end else "out={0}".format(*output)
    log:
        "{sample}/logs/{sample}_quality_filter.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    resources:
        mem = config.get("java_mem", JAVA_MEM)
    shell:
        """{SHPFXM} bbduk2.sh {params.inputs} \
              {params.outputs} \
               {params.rref} {params.lref} \
               {params.mink} qout=33 stats={output.stats} \
               {params.hdist} {params.k} trimq={params.trimq} \
               qtrim={params.qtrim} threads={threads} \
               minlength={params.minlength} trd=t \
               minbasefrequency={params.minbasefrequency} \
               interleaved={params.interleaved}\
               overwrite=true \
               ecco={params.error_correction_pe} \
               -Xmx{resources.mem}G 2> {log}
        """
# if there are no references, decontamination will be skipped

if len(config.get("contaminant_references", {}).keys()) > 0:

    processed_steps += ['clean']

    rule build_decontamination_db:
        output:
            "ref/genome/1/summary.txt"
        threads:
            config.get("threads", 1)
        resources:
            mem = config.get("java_mem", JAVA_MEM)
        log:
            "logs/build_decontamination_db.log"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        params:
            k = config.get("contaminant_kmer_length", CONTAMINANT_KMER_LENGTH),
            refs_in = " ".join(["ref_%s=%s" % (n, fa) for n,fa in config["contaminant_references"].items()]),
        shell:
            """{SHPFXM} bbsplit.sh -Xmx{resources.mem}G {params.refs_in} threads={threads} k={params.k} local=t 2> {log}"""


    rule decontamination:
        input:
            expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                step=processed_steps[-2], fraction=multifile_fractions),
            db = "ref/genome/1/summary.txt"
        output:
            temp(expand("{{sample}}/sequence_quality_control/{{sample}}_{step}_{fraction}.fastq.gz",
                fraction=multifile_fractions, step=processed_steps[-1])),
            contaminants = expand("{{sample}}/sequence_quality_control/contaminants/{db}_{fraction}.fastq.gz",
                    db=list(config["contaminant_references"].keys()),
                    fraction=multifile_fractions),
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
            paired= "true" if paired_end else "false",
            input_single= lambda wc,input: input[2] if paired_end else input[0],
            output_single= lambda wc,output: output[2] if paired_end else output[0]
        log:
            "{sample}/logs/{sample}_decontamination.log"
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads:
            config.get("threads", 1)
        resources:
            mem = config.get("java_mem", JAVA_MEM)
        shell:
            """
            if [ "{params.paired}" = true ] ; then
            {SHPFXM} bbsplit.sh in1={input[0]} in2={input[1]} \
                    outu1={output[0]} outu2={output[1]} \
                    basename="{params.contaminant_folder}/%_R#.fastq.gz" \
                    maxindel={params.maxindel} minratio={params.minratio} \
                    minhits={params.minhits} ambiguous={params.ambiguous} refstats={output.stats}\
                    threads={threads} k={params.k} local=t \
                    -Xmx{resources.mem}G 2> {log}
            fi

            {SHPFXM} bbsplit.sh in={params.input_single}  \
                    outu={params.output_single} \
                   basename="{params.contaminant_folder}/%_se.fastq.gz" \
                   maxindel={params.maxindel} minratio={params.minratio} \
                   minhits={params.minhits} ambiguous={params.ambiguous} refstats={output.stats} append \
                   interleaved=f threads={threads} k={params.k} local=t \
                   -Xmx{resources.mem}G 2>> {log}


            """



processed_steps += ['QC']

def get_ribosomal_rna_input(wildcards):

    inputs = []
    data_type = config["samples"][wildcards.sample].get("type", "metagenome").lower()

    clean_reads = "{sample}/sequence_quality_control/{sample}_{step}_{fraction}.fastq.gz".format(step=processed_steps[-2],**wildcards)
    rrna_reads = "{sample}/sequence_quality_control/contaminants/rRNA_{fraction}.fastq.gz".format(**wildcards)

    if data_type == "metagenome" and os.path.exists(rrna_reads):
        return [clean_reads, rrna_reads]
    else:
        return [clean_reads]


rule postprocess_after_decontamination:
    input:
        get_ribosomal_rna_input
    output:
        "{{sample}}/sequence_quality_control/{{sample}}_{step}_{{fraction}}.fastq.gz".format(step=processed_steps[-1])
    threads:
        1
    shell:
        "{SHPFXS} cat {input} > {output}"




if paired_end:
    rule calculate_insert_size:
        input:
            unpack(get_quality_controlled_reads)
        output:
            ihist = "{sample}/sequence_quality_control/read_stats/QC_insert_size_hist.txt",
            cardinality= temp("{sample}/sequence_quality_control/read_stats/QC_cardinality.txt"),
            read_length= "{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt"
        threads:
            config.get("threads", 1)
        resources:
            mem = config.get("java_mem", JAVA_MEM)
        conda:
            "%s/required_packages.yaml" % CONDAENV
        log:
            "{sample}/logs/{sample}_calculate_insert_size.log"
        benchmark:
            "logs/benchmarks/merge_pairs/{sample}_insert_size.txt"
        params:
            kmer = config.get("merging_k", MERGING_K),
            extend2 = config.get("merging_extend2", MERGING_EXTEND2),
            flags = 'loose ecct'
        shell:
            """bbmerge.sh -Xmx{resources.mem}G threads={threads} \
                   in1={input.R1} in2={input.R2} \
                   {params.flags} k={params.kmer} \
                   extend2={params.extend2} \
                   ihist={output.ihist} outc={output.cardinality} merge=f \
                   mininsert0=35 minoverlap0=8 2> >(tee {log})
                
                readlength.sh in={input.R1} in2={input.R2} out={output.read_length} 2> >(tee {log})
            """
else:

    rule calculate_read_length_hist:
        input:
            unpack(get_quality_controlled_reads)
        output:
            read_length= "{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt",
            cardinality= temp("{sample}/sequence_quality_control/read_stats/QC_cardinality.txt")
        params:
            kmer = config.get("merging_k", MERGING_K),
        threads:
            config.get("threads", 1)
        resources:
            mem = config.get("java_mem", JAVA_MEM)
        conda:
            "%s/required_packages.yaml" % CONDAENV
        log:
            "{sample}/logs/{sample}_calculate_read_length.log"
        shell:
            """ 
                readlength.sh in={input.se} out={output.read_length} 2> >(tee {log})
                loglog.sh in={input.se} k={params.kmer} > {output.cardinality} 2> >(tee {log})
            """

localrules: combine_read_length_stats, combine_insert_stats, combine_cardinality, combine_read_counts


rule combine_read_length_stats:
    input:
        expand("{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt",sample=SAMPLES),
    output:
        'stats/read_length_stats.tsv'
    run:
        import pandas as pd
        import os



        Stats= pd.DataFrame()

        for length_file in input:
            sample= length_file.split(os.path.sep)[0]

            data = parse_comments(length_file)
            data = pd.Series(data)[['Reads','Bases','Max','Min','Avg','Median','Mode','Std_Dev']]

            Stats[sample]=data

        Stats.to_csv(output[0],sep='\t')   

             
rule combine_cardinality:
    input:
        expand("{sample}/sequence_quality_control/read_stats/QC_cardinality.txt",sample=SAMPLES),
    output:
        'stats/cardinality.tsv'
    run:
        import pandas as pd
        import os

        Stats= pd.Series()

        for file in input:
            sample= file.split(os.path.sep)[0]
            with open(file) as f:
                cardinality= int(f.read().strip())

            Stats.loc[sample]=cardinality

        Stats.to_csv(output[0],sep='\t')



if paired_end:
    rule combine_insert_stats:
        input:
            expand("{sample}/sequence_quality_control/read_stats/QC_insert_size_hist.txt",sample=SAMPLES),
        output:
            'stats/insert_stats.tsv'
        run:
            import pandas as pd
            import os
            Stats= pd.DataFrame()

            for insert_file in input:
                sample= insert_file.split(os.path.sep)[0]

                data = parse_comments(insert_file)
                data = pd.Series(data)[['Avg','Median','Mode','STDev','PercentOfPairs']]

                Stats[sample]=data

            Stats.T.to_csv(output[0],sep='\t')  


rule combine_read_counts:
    input:
        expand("{sample}/sequence_quality_control/read_stats/read_counts.tsv",sample=SAMPLES)
    output:
        "stats/read_counts.tsv"
    run:
        import pandas as pd

        Read_stats=pd.DataFrame()

        for f in input:
            d= pd.read_table(f,index_col=[0,1])
            Read_stats=Read_stats.append(d)

        Read_stats.T.to_csv(output[0],sep='\t')


rule finalize_QC:
    input:
        unpack(get_quality_controlled_reads),
            rules.decontamination.output.contaminants,
            "{sample}/sequence_quality_control/{sample}_decontamination_reference_stats.txt",
            "{sample}/logs/{sample}_quality_filtering_stats.txt",
            expand("{{sample}}/sequence_quality_control/read_stats/{step}.zip", step=processed_steps),
            read_count_files= expand("{{sample}}/sequence_quality_control/read_stats/{step}_read_counts.tsv", step=processed_steps),
            read_length_hist="{sample}/sequence_quality_control/read_stats/QC_read_length_hist.txt"

    output:
        touch("{sample}/sequence_quality_control/finished_QC"),
        read_stats= "{sample}/sequence_quality_control/read_stats/read_counts.tsv" # exists alredy before
    run:
        print("Finished QC for sample {sample}\n".format(**wildcards))
        import pandas as pd
        All_read_counts= pd.DataFrame()
        for read_stats_file in input.read_count_files:
            d= pd.read_table(read_stats_file,index_col=[0,1])
            All_read_counts= All_read_counts.append(d)
        All_read_counts.to_csv(output.read_stats,sep='\t')


# 
rule QC_report:
    input:
        expand("{sample}/sequence_quality_control/finished_QC",sample=SAMPLES),
        "stats/read_counts.tsv",
        'stats/cardinality.tsv',
        read_length_stats= ['stats/insert_stats.tsv','stats/read_length_stats.tsv'] if paired_end else 'stats/read_length_stats.tsv'
    output:
        touch("finished_QC")
    shell:
        """
        if [ -d ref ]; then
            rm -r ref
        fi
        """


    # aggregate stats reports ...
