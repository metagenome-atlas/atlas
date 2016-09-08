import os
from subprocess import check_output


def read_count(fastq):
    total = 0
    count_file = fastq + '.count'
    if os.path.exists(fastq) and os.path.getsize(fastq) > 100:
        if not os.path.exists(count_file):
            check_output("awk '{n++}END{print n/4}' %s > %s" % (fastq, fastq + '.count'), shell=True)
        with open(count_file) as fh:
            for line in fh:
                total = int(line.strip())
                break
    return total


def get_samples(eid, dir="demultiplexed", coverage_cutoff=1000):
    """Grab samples from files residing in results/<eid>/<dir>. Expecting files with .fastq and a
    naming convention like <sample>_R1.fastq and <sample>_R2.fastq.
    """
    samples = set()
    input_dir = os.path.join("results", eid, "demux")
    for f in os.listdir(input_dir):
        if f.endswith("fastq") and ("_r1" in f or "_R1" in f):
            if read_count(os.path.join(input_dir, f)) > coverage_cutoff:
                samples.add(f.partition(".")[0].partition("_")[0])
    return samples


EID = config['eid']
SAMPLES = get_samples(EID)
CONTAMINANT_REFS = config['contamination_filtering']['references'].split(",")


rule all:
    input:
        # desired output files to keep


rule build_contaminant_references:
    # check for file compression
    input: "ref/contamination_references/{fasta}"
    output:
        f1 = "{fasta}.1.bt2",
        f2 = "{fasta}.2.bt2",
        f3 = "{fasta}.3.bt2",
        f4 = "{fasta}.4.bt2",
        r1 = "{fasta}.rev.1.bt2",
        r2 = "{fasta}.rev.2.bt2"
    shell: "bowtie2-build {input} {input}"


rule quality_filter_reads:
    input:
        r1 = "results/{eid}/demultiplexed/{sample}_R1.fastq",
        r2 = "results/{eid}/demultiplexed/{sample}_R2.fastq"
    output:
        r1 = "results/{eid}/demultiplexed/filtered/{sample}_filtered_R1.fastq",
        r2 = "results/{eid}/demultiplexed/filtered/{sample}_filtered_R2.fastq",
        stats = temp("results/{eid}/demultiplexed/filtered/{sample}_quality_filtering_stats.txt")
    message: "Filtering reads using BBDuk2 to remove adapters and phiX with matching kmer length of {params.k} at a hamming distance of {params.hdist} and quality trim both ends to Q{params.quality}. Reads shorter than {params.minlength} were discarded."
    params:
        adapters = config['filtering']['adapters'],
        quality = config['filtering']['minimum_base_quality'],
        hdist = config['filtering']['allowable_kmer_mismatches'],
        k = config['filtering']['reference_kmer_match_length'],
        qtrim = "rl",
        ktrim = "l",
        minlength = config['filtering']['minimum_passing_read_length']
    threads: config['filtering']['threads']
    shell: """bbduk2.sh -Xmx8g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} \
               fref={params.adapters} stats={output.stats} hdist={params.hdist} k={params.k} \
               trimq={params.quality} qtrim={params.qtrim} threads={threads} ktrim={params.ktrim} \
               minlength={params.minlength} overwrite=true
           """


rule count_filtered_reads:
    input: rules.quality_filter_reads.output.r1
    output: "results/{eid}/logs/{sample}_filtered_R1.fastq.count"
    threads: 1
    shell: "awk '{{n++}}END{{print n/4}}' {input} > {output}"


rule join_reads:
    input:
        r1 = rules.quality_filter_reads.r1,
        r2 = rules.quality_filter_reads.r2
    output:
        joined = "results/{eid}/joined/{sample}.extendedFrags.fastq",
        hist = "results/{eid}/joined/{sample}.hist",
        failed_r1 = "results/{eid}/joined/{sample}.notCombined_1.fastq",
        failed_r2 = "results/{eid}/joined/{sample}.notCombined_2.fastq"
    message: "Joining reads using `flash`"
    shadow: "shallow"
    params:
        min_overlap = config['merging']['minimum_overlap'],
        max_overlap = config['merging']['maximum_overlap'],
        max_mismatch_density = config['merging']['maximum_mismatch_density'],
        phred_offset = config['phred_offset']
    log: "results/{eid}/logs/{sample}_flash.log"
    threads: config['merging']['flash_threads']
    shell: """flash {input.r1} {input.r2} --min-overlap {params.min_overlap} \
                  --max-overlap {params.max_overlap} --max-mismatch-density {params.max_mismatch_density} \
                  --phred-offset {params.phred_offset} --output-prefix {} \
                  --output-directory {} --threads {threads}
           """
