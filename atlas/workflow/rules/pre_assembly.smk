
def get_preprocessing_steps(config):
    preprocessing_steps = ["QC"]
    if config.get("normalize_reads_before_assembly", False):
        preprocessing_steps.append("normalized")

    if config.get("error_correction_before_assembly", True):
        preprocessing_steps.append("errorcorr")

    if config.get("merge_pairs_before_assembly", True) and PAIRED_END:
        preprocessing_steps.append("merged")

    return ".".join(preprocessing_steps)


assembly_preprocessing_steps = get_preprocessing_steps(config)

# I have problems with se reads
if SKIP_QC & (len(MULTIFILE_FRACTIONS) < 3):

    rule init_pre_assembly_processing:
        input:  #expect SE or R1,R2 or R1,R2,SE
            get_quality_controlled_reads,
        output:
            temp(
                expand(
                    "{{sample}}/assembly/reads/QC_{fraction}.fastq.gz",
                    fraction=MULTIFILE_FRACTIONS,
                )
            ),
        params:
            inputs=lambda wc, input: io_params_for_tadpole(input, "in"),
            interleaved=(
                lambda wc: "t"
                if (config.get("interleaved_fastqs", False) & SKIP_QC)
                else "f"
            ),
            outputs=lambda wc, output: io_params_for_tadpole(output, "out"),
            verifypaired="t" if PAIRED_END else "f",
        log:
            "{sample}/logs/assembly/init.log",
        conda:
            "%s/required_packages.yaml" % CONDAENV
        threads: config.get("simplejob_threads", 1)
        resources:
            mem=config["simplejob_mem"],
            java_mem=int(config["simplejob_mem"] * JAVA_MEM_FRACTION),
        shell:
            """
            reformat.sh {params.inputs} \
                interleaved={params.interleaved} \
                {params.outputs} \
                iupacToN=t \
                touppercase=t \
                qout=33 \
                overwrite=true \
                verifypaired={params.verifypaired} \
                addslash=t \
                trimreaddescription=t \
                threads={threads} \
                pigz=t unpigz=t \
                -Xmx{resources.java_mem}G 2> {log}
            """


else:

    localrules:
        init_pre_assembly_processing,

    rule init_pre_assembly_processing:
        input:
            lambda wildcards: get_quality_controlled_reads(wildcards, include_se=True),
        output:
            temp(
                expand(
                    "{{sample}}/assembly/reads/QC_{fraction}.fastq.gz",
                    fraction=MULTIFILE_FRACTIONS,
                )
            ),
        threads: 1
        run:
            # make symlink
            assert len(input) == len(
                output
            ), "Input and ouput files have not same number, can not create symlinks for all."
            for i in range(len(input)):
                os.symlink(os.path.abspath(input[i]), output[i])




#
rule normalize_reads:
    input:
        expand(
            "{{sample}}/assembly/reads/{{previous_steps}}_{fraction}.fastq.gz",
            fraction=MULTIFILE_FRACTIONS,
        ),
    output:
        reads=temp(
            expand(
                "{{sample}}/assembly/reads/{{previous_steps}}.normalized_{fraction}.fastq.gz",
                fraction=MULTIFILE_FRACTIONS,
            )
        ),
        histin="{sample}/assembly/normalization/histogram_{previous_steps}_before_normalization.tsv.gz",
        histout=(
            "{sample}/assembly/normalization/histogram_{previous_steps}_after.tsv.gz"
        ),
    params:
        k=config.get("normalization_kmer_length", NORMALIZATION_KMER_LENGTH),
        target=config.get("normalization_target_depth", NORMALIZATION_TARGET_DEPTH),
        mindepth=config["normalization_minimum_kmer_depth"],
        inputs=lambda wc, input: io_params_for_tadpole(input),
        outputs=lambda wc, output: io_params_for_tadpole(output.reads, key="out"),
        tmpdir="tmpdir=%s" % TMPDIR if TMPDIR else "",
    log:
        "{sample}/logs/assembly/pre_process/normalization_{previous_steps}.log",
    benchmark:
        "logs/benchmarks/assembly/pre_process/normalization/{sample}_{previous_steps}.txt"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    threads: config.get("threads", 1)
    resources:
        mem=config["mem"],
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
    shell:
        " bbnorm.sh {params.inputs} "
        " {params.outputs} "
        " {params.tmpdir} "
        " tossbadreads=t "
        " hist={output.histin} "
        " histout={output.histout} "
        " mindepth={params.mindepth} "
        " k={params.k} "
        " target={params.target} "
        " prefilter=t "
        " threads={threads} "
        " -Xmx{resources.java_mem}G &> {log} "


rule error_correction:
    input:
        expand(
            "{{sample}}/assembly/reads/{{previous_steps}}_{fraction}.fastq.gz",
            fraction=MULTIFILE_FRACTIONS,
        ),
    output:
        temp(
            expand(
                "{{sample}}/assembly/reads/{{previous_steps}}.errorcorr_{fraction}.fastq.gz",
                fraction=MULTIFILE_FRACTIONS,
            )
        ),
    benchmark:
        "logs/benchmarks/assembly/pre_process/{sample}_error_correction_{previous_steps}.txt"
    log:
        "{sample}/logs/assembly/pre_process/error_correction_{previous_steps}.log",
    conda:
        "%s/required_packages.yaml" % CONDAENV
    resources:
        mem=config["mem"],
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
    params:
        inputs=lambda wc, input: io_params_for_tadpole(input),
        outputs=lambda wc, output: io_params_for_tadpole(output, key="out"),
        prefilter=2,  # Ignore kmers with less than 2 occurance
        minprob=config["error_correction_minprob"],
        tossdepth=config["error_correction_minimum_kmer_depth"],
        tossjunk="t" if config["error_correction_remove_lowdepth"] else "f",
        lowdepthfraction=config["error_correction_lowdepth_fraction"],
        aggressive=config["error_correction_aggressive"],
        shave="f",  # Shave and rinse can produce substantially better assemblies for low-depth data, but they are very slow for large metagenomes.
    threads: config.get("threads", 1)
    shell:
        "tadpole.sh -Xmx{resources.java_mem}G "
        " prefilter={params.prefilter} "
        " prealloc=1 "
        " {params.inputs} "
        " {params.outputs} "
        " mode=correct "
        " aggressive={params.aggressive} "
        " tossjunk={params.tossjunk} "
        " lowdepthfraction={params.lowdepthfraction}"
        " tossdepth={params.tossdepth} "
        " merge=t "
        " shave={params.shave} rinse={params.shave} "
        " threads={threads} "
        " pigz=t unpigz=t "
        " ecc=t ecco=t "
        "&> {log} "


rule merge_pairs:
    input:
        expand(
            "{{sample}}/assembly/reads/{{previous_steps}}_{fraction}.fastq.gz",
            fraction=["R1", "R2"],
        ),
    output:
        temp(
            expand(
                "{{sample}}/assembly/reads/{{previous_steps}}.merged_{fraction}.fastq.gz",
                fraction=["R1", "R2", "me"],
            )
        ),
    threads: config.get("threads", 1)
    resources:
        mem=config["mem"],
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
    conda:
        "%s/required_packages.yaml" % CONDAENV
    log:
        "{sample}/logs/assembly/pre_process/merge_pairs_{previous_steps}.log",
    benchmark:
        "logs/benchmarks/assembly/pre_process/merge_pairs_{previous_steps}/{sample}.txt"
    shadow:
        "shallow"
    params:
        kmer=config.get("merging_k", MERGING_K),
        extend2=config.get("merging_extend2", MERGING_EXTEND2),
        flags=config.get("merging_flags", MERGING_FLAGS),
    shell:
        """
        bbmerge.sh -Xmx{resources.java_mem}G threads={threads} \
            in1={input[0]} in2={input[1]} \
            outmerged={output[2]} \
            outu={output[0]} outu2={output[1]} \
            {params.flags} k={params.kmer} \
            pigz=t unpigz=t \
            extend2={params.extend2} 2> {log}
        """
