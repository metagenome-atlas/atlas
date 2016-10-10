import os
import glob
from subprocess import check_output
from util.IO import cat_reads
from util.IO import interleave_reads


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
    input_dir = os.path.join("results", eid, dir)
    for f in os.listdir(input_dir):
        if f.endswith("fastq") and ("_r1" in f or "_R1" in f):  # could probably replex with a regex
            if read_count(os.path.join(input_dir, f)) > coverage_cutoff:
                samples.add(f.partition(".")[0].partition("_")[0])
    return samples


def get_databases():
    """Grab databases from files residing in /databases/taxonomic/, /databases/functional/,
    and /databases/contaminant/. Expecting a directory where the only files are database files."""

    # gather files
    tax = glob.glob('./databases/taxonomic/*')
    func = glob.glob('./databases/functional/*')
    contam = glob.glob('./databases/contaminant/*')

    # filter non files
    tax = [os.path.basename(x) for x in tax if os.path.isfile(x)]
    func = [os.path.basename(x) for x in func if os.path.isfile(x)]
    contam = [os.path.basename(x) for x in contam if os.path.isfile(x)]

    # compose dictionary
    dbs = {'taxonomic': tax, 'functional': func, 'contaminants': contam}

    return dbs



EID = config['eid']
SAMPLES = get_samples(EID)
DBS = get_databases():


rule all:
    input:
        # desired output files to keep


rule gunzip:
    input:
        zipfile = "input/{eid}/{sample}."


rule build_contaminant_references:
    input:
        contaminant_db = "databases/contaminant/{db}",

        # can I call Python code here?
        contaminant_db_name = os.path.splitext(contaminant_db)
    output:
        f1 = "{fasta}.1.bt2",
        f2 = "{fasta}.2.bt2",
        f3 = "{fasta}.3.bt2",
        f4 = "{fasta}.4.bt2",
        r1 = "{fasta}.rev.1.bt2",
        r2 = "{fasta}.rev.2.bt2"
    message:
        "Formatting contaminant databases"
    shell:
        "bowtie2-build {input.contaminant_db} {input.contaminant_db_name}"


rule build_functional_databases:
    input:
        functional_db = "databases/functional/{db}"
    output:
        f1 = "{lastal_database}.bck",
        f2 = "{lastal_database}.des",
        f3 = "{lastal_database}.prj",
        f4 = "{lastal_database}.sds",
        f5 = "{lastal_database}.ssp",
        f6 = "{lastal_database}.suf",
        f7 = "{lastal_database}.tis",
        f8 = "{lastal_database}-names.txt"
    message:
        "Formatting functional databases"
    shell:
        "lastdb+ {input.taxonomic_db} {input.taxonomic_db} -p"


rule build_taxonomic_databases:
    input:
        taxonomic_db = "annotation_dbs/taxonomic/{db}"
    output:
        f1 = "{lastal_database}.bck",
        f2 = "{lastal_database}.des",
        f3 = "{lastal_database}.prj",
        f4 = "{lastal_database}.sds",
        f5 = "{lastal_database}.ssp",
        f6 = "{lastal_database}.suf",
        f7 = "{lastal_database}.tis",
        f8 = "{lastal_database}-names.txt"
    message:
        "Formatting taxonomic databases"
    shell:
        "lastdb+ {input.taxonomic_db} {input.taxonomic_db}"


rule join_reads:
    input:
        r1 = rules.quality_filter_reads.output.r1,
        r2 = rules.quality_filter_reads.output.r2
    output:
        joined = "results/{eid}/joined/{sample}.extendedFrags.fastq",
        hist = "results/{eid}/joined/{sample}.hist",
        failed_r1 = "results/{eid}/joined/{sample}.notCombined_1.fastq",
        failed_r2 = "results/{eid}/joined/{sample}.notCombined_2.fastq"
    message:
        "Joining reads using `flash`"
    shadow:
        "shallow"
    params:
        min_overlap = config['merging']['minimum_overlap'],
        max_overlap = config['merging']['maximum_overlap'],
        max_mismatch_density = config['merging']['maximum_mismatch_density'],
        phred_offset = config['phred_offset']
    log:
        "results/{eid}/logs/{sample}_flash.log"
    threads:
        config['merging']['flash_threads']
    shell:
        """flash {input.r1} {input.r2} --min-overlap {params.min_overlap} \
                  --max-overlap {params.max_overlap} --max-mismatch-density {params.max_mismatch_density} \
                  --phred-offset {params.phred_offset} --output-prefix {wildcards.sample} \
                  --output-directory results/{wildcards.eid}/joined/ --threads {threads}
        """


rule filter_contaminants:
    input:
        joined = rules.join_reads.output.joined,
        r1 = rules.join_reads.output.failed_r1,
        r2 = rules.join_reads.output.failed_r2,
        prefix = rules.combine_contaminant_references.output
    output:
        "results/{eid}/joined/{sample}_joined_filtered.fastq"
    message:
        "Aligning all joined and reads that failed to join as single-end against the contamination reference."
    shadow:
        "shallow"
    threads:
        config['contamination_filtering']['threads']
    shell:
        """bowtie2 --threads {threads} --very-sensitive-local -x {input.prefix} -q -U {input.joined},{input.r1},{input.r2} \
                  | samtools view -@ {threads} -hf4 \
                  | samtools sort -@ {threads} -T {wildcards.sample} -o -m 8G - \
                  | bedtools bamtofastq -i stdin -fq {output}
           """


rule trim_reads:
    input:
        filtered = rule.filter_contaminants.output,
        picard = ?
    output:
        joined = "results/{eid}/trimmed/{sample}_trimmed_filtered_joined.fastq",
        R1 = "results/{eid}/trimmed/{sample}_trimmed_filtered_R1.fastq",
        R2 = "results/{eid}/trimmed/{sample}_trimmed_filtered_R2.fastq"
    params:
        single_end = config['SE'],
        phred_value = config['phred33'],
        min_length = config['length?']
    message:
        "Trimming filtered reads using trimmomatic"
    log:
        "results/{eid}/trimmed/{sample}.log"
    shell:
        """java -Xmx32g -jar trimmomatic-0.33.jar SE -phred33 {input.filtered} -trimlog {log} \
            ILLUMINACLIP:adapters/TruSeq2-SE:2:30:10 LEADING:3 TRAILING:3 \
            SLIDINGWINDOW:4:15 MINLEN:{params.min_length} {output}
        """


rule fastqc_R1:
    input:
        rule.trim_reads.output.R1
    output:
        "results/{eid}/qc/{sample}_R1.fastq"
    shell:
        "fastqc {input} -o {output}"


rule fastqc_R2:
    input:
        rule.trim_reads.output.R2
    output:
        "results/{eid}/qc/{sample}_R2.fastq"
    shell:
        "fastqc {input} -o {output}"


rule fastqc_joined:
    input:
        rule.trim_reads.output.joined
    output:
        "results/{eid}/qc/{sample}_joined.fastq"
    shell:
        "fastqc {input} -o {output}"


# For metatranscriptomes only!
rule remove_rRNAs:
    input:
        rule.interleave_reads.output,
        rule.trim_reads.joined.output
    output:
        joined = "results/{eid}/trimmed/{sample}_trimmed_filtered_no-rrna_joined.fastq",
        R1 = "results/{eid}/trimmed/{sample}_trimmed_filtered_no-rrna_R1.fastq",
        R2 = "results/{eid}/trimmed/{sample}_trimmed_filtered_no-rrna_R2.fastq"
    message:
        "Aligning all joined and reads to remove rRNAs for mRNA de novo transcriptome assembly."
    shadow:
        "shallow"
    threads:
        config['contamination_filtering']['threads']
    shell:
        """bowtie2 --threads {threads} --very-sensitive-local -x {input.prefix} -q -U {input.joined},{input.r1},{input.r2} \
              | samtools view -@ {threads} -hf4 \
              | samtools sort -@ {threads} -T {wildcards.sample} -o -m 8G - \
              | bedtools bamtofastq -i stdin -fq {output}
        """


rule interleave_reads:
    input:
        r1 = rules.join_reads.output.failed_r1,
        r2 = rules.join_reads.output.failed_r2
    output:
        "results/{eid}/trimmed/{sample}_trimmed_filtered_interleaved.fastq"
    message:
        "Interleaving non combined R1 and R2 reads"
    run:
        IO.interleave_reads(input.r1, input.r2, output)


rule merge_joined_interleaved:
    input:
        il = rules.interleave_reads.output,
        joined = rules.join_reads.output.joined
    output:
        "results/{eid}/assembly/{sample}_merged.fastq"
    message:
        "Merging joined and interleaved reads"
    run:
        IO.cat_reads(input.il, input.joined, output)


rule annotate_reads_rRNA:
    input:
        rule.interleave_reads.output,
        trim_reads = rule.trim_reads.joined.output
    output:
        "results/{eid}/annotation/reads/{sample}_{database}"
    message:
        "Annotation of rRNAs in quality controlled reads"
    params:
        top_hit = config['lastplus']['top_best_hit'],
        e_value_cutoff = config['lastplus']['e_value_cutoff'],
        bit_score_cutoff = config['lastplus']['bit_score_cutoff']
    threads:
        config['lastplus']['threads']
    shell:
        """lastal+ -P {threads} -K {params.top_hit} -E {params.e_value_cutoff} -S {params.bit_score_cutoff} -o {output} \
            {input.database} {input.trim_reads}"""


rule taxonomic_placement_rRNAs_LCA:
    input:
        rule.annotate_reads_rRNA.output
    output:
        "results/{eid}/annotation/reads/{sample}_{database}"
    message:
        "Parse rRNA reads and place taxonomy using LCA++"
    params:
        # TODO
    threads:
        config['lastplus']['threads']
    shell:


# will want to change this as we add assemblers
rule assemble:
    input:
        # TODO
    output:
        # TODO


# for metagenomes only
rule megahit:
    input:
        rules.filter_contaminants.output
    output:
        "results/{eid}/assembly/{sample}.contigs.fa"
    params:
        memory = config['assembly']['memory'],
        min_count = config['assembly']['minimum_count'],
        k_min = config['assembly']['kmer_min'],
        k_max = config['assembly']['kmer_max'],
        k_step = config['assembly']['kmer_step'],
        merge_level = config['assembly']['merge_level'],
        prune_level = config['assembly']['prune_level'],
        low_local_ratio = config['assembly']['low_local_ratio'],
        min_contig_len = config['assembly']['minimum_contig_length']
    message:
        "Assembling using megahit"
    threads:
        config['assembly']['threads']
    shell:
        """megahit --num-cpu-threads {threads} --memory {params.memory} --read {input} \
        --k-min {params.k_min} --k-max {params.k_max} --k-step {params.k_step} \
        --out-dir results/{wildcards.eid}/assembly --out-prefix {wildcards.sample} \
        --min-contig-len {params.min_contig_len} --min-count {params.min_count} \
        --merge-level {params.merge_level} --prune-level {params.prune_level} \
        --low-local-ratio {params.low_local_ratio}"""


# for metatranscriptomes only
rule trinity:
    input:
        # need to replace with reference to rule
        extendedFrags = 'results/{eid}/trimmed/{sample}_trimmed_extendedFrags.fastq',  # after we fix trimming!
        interleaved = 'results/{eid}/interleaved/{sample}_trimmed_interleaved.fastq'
    output:
        "results/{eid}/assembly/{sample}.contigs.fa"
    params:
        seqtype = config['assembly']['seqtype'],  # default fastq
        read_pairing = config['assembly']['single'],  # default single for extendedFrags
        memory = config['assembly']['max_memory'],
        run_as_paired = config['assembly']['run_as_paired']
    message:
        "Assembling using trinity"
    threads:
        config['assembly']['threads']
    shell:
        """Trinity --seqType {params.seqtype} --single {input.extendedFrags}, {input.interleaved} \
        --run_as_paired --max_memory {params.max_memory} --CPU {threads}"""


rule length_filter:
    input:
        rules.assemble.output
    output:
        passing = "results/{eid}/assembly/{sample}_length_pass.fa",
        fail = "results/{eid}/assembly/{sample}_length_fail.fa"
    params:
        min_contig_length = config['assembly']['filtered_contig_length']
    shell:
        """python scripts/fastx.py length-filter --min-length {params.min_contig_length} \
        {input} {output.passing} {output.fail}"""


rule assembly_stats:
    input:
        assembled = rules.assembly.output,
        filtered = rules.length_filter.output
    output:
        assembled = "results/{eid}/assembly/{sample}_length_fail_assembly-stats.txt",
        filtered = "results/{eid}/assembly/{sample}_length_pass_assembly-stats.txt"
    message:
        "Obtaining assembly statistics"
    shell:
        """perl scripts/CountFasta.pl {input.assembled} > {output.assembled}
           perl scripts/CountFasta.pl {input.filtered} > {output.filtered}
        """


# have Joe review
rule fgsplus_passed:
    input:
        rules.length_filter.output.passing,
    output:
        "results/{eid}/annotation/orfs/{sample}_length_pass.faa",
    params:
        sem = config['annotation']['sequencing_error_model'],
        memory = config['annotation']['memory']
    threads:
        config['annotation']['threads']
    shell:
        "FGS+ -s {input} -o {output} -w 1 -t {params.sem} -p {threads} -m {params.memory}"


rule fgsplus_failed:
    input:
        rules.length_filter.output.fail
    output:
        "results/{eid}/annotation/orfs/{sample}_length_fail.faa"
    params:
        sem = config['annotation']['sequencing_error_model'],
        memory = config['annotation']['memory']
    threads:
        config['annotation']['threads']
    shell:
        "FGS+ -s {input} -o {output} -w 1 -t {params.sem} -p {threads} -m {params.memory}"


rule prodigal_orfs:
    input:
        rules.length_filter.output.passing
    output:
        prot = "results/{eid}/annotation/orfs/{sample}.faa",
        nuc = "results/{eid}/annotation/orfs/{sample}.fasta",
        gff = "results/{eid}/annotation/orfs/{sample}.gff"
    params:
        g = config['annotation']['translation_table']
    shell:
        "prodigal -i {input} -o {output.gff} -f gff -a {output.prot} -d {output.nuc} -g {params.g} -p meta"


rule maxbin_bins:
    input:
        reads = rules.filter_contaminants.output
        contigs = rules.assemble.output
    output:
        bins = "results/{eid}/binning/{sample}.fasta",
        abundance = "results/{eid}/binning/{sample}.abund1",
        log = "results/{eid}/binning/{sample}.log",
        marker = "results/{eid}/binning/{sample}.marker",
        summary = "results/{eid}/binning/{sample}.summary",
        tooshort = "results/{eid}/binning/{sample}.tooshort",
        noclass = "results/{eid}/binning/{sample}.noclass"
    params:
        min_contig_len = config['binning']['minimum_contig_length'],
        max_iteration = config['binning']['maximum_iterations'],
        prob_threshold = config['binning']['probability_threshold'],
        markerset = config['binning']['marker_set']
    threads:
        config['binning']['threads']
    shell:
        """run_MaxBin.pl -contig {input.contigs} -out {output} -reads {input.reads} \
        -min_contig_length {params.min_contig_len} -max_iteration {params.max_iteration} \
        -thread {threads} -markerset {params.markerset}"""


rule lastplus_orfs
    input:  # how to do wrap rule or ifelse?
        fgsplus_orfs = rules.fgsplus.output,
        prodigal_orfs = rules.prodigal.output,
        database = rules.format_database.output
    output:
        annotation = "results/{eid}/annotation/orfs/{sample}_{database}"
    params:
        top_hit = config['lastplus']['top_best_hit'],
        e_value_cutoff = config['lastplus']['e_value_cutoff'],
        bit_score_cutoff = config['lastplus']['bit_score_cutoff']
    threads:
        config['lastplus']['threads']
    shell:
        """lastal+ -P {threads} -K {params.top_hit} -E {params.e_value_cutoff} -S {params.bit_score_cutoff} -o {output} \
        {input.database} {input.fgsplus_orfs}"""
