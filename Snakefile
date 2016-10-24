import os
from glob import glob
from subprocess import check_output
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


def get_samples(path, coverage_cutoff=1000):
    """Expecting files with .fastq and a naming convention like <sample>_R1.fastq and
    <sample>_R2.fastq.
    """
    samples = set()
    for f in os.listdir(path):
        if f.endswith("fastq") and ("_r1" in f or "_R1" in f):
            if read_count(os.path.join(path, f)) > coverage_cutoff:
                samples.add(f.partition(".")[0].partition("_")[0])
    return samples


def pattern_search(path, patterns):
    """Grab files under path that match pattern."""
    files = []
    for pattern in patterns:
        files.extend(glob("%s/%s" % (path, pattern)))
    return [os.path.basename(i) for i in files]


# snakemake --configfile config/atlas_config.yaml --config eid=test-experiment
configfile: "./config/atlas_config.yaml"

EID = config['eid']
SAMPLES = get_samples(os.path.join("input", EID), 200)
CONTAMINANT_DBS = pattern_search("databases/contaminant", ["*.fa", "*.fasta"])
FUNCTIONAL_DBS = pattern_search("databases/functional", ["*.fa", "*.fasta"])
TAXONOMIC_DBS = pattern_search("databases/taxonomic", ["*.fa", "*.fasta"])


rule all:
    input:
        # these can eventually be removed. testing purposes only.
        # contaminants
        expand("databases/contaminant/{db}.1.bt2", db=CONTAMINANT_DBS),
        expand("databases/contaminant/{db}.2.bt2", db=CONTAMINANT_DBS),
        expand("databases/contaminant/{db}.3.bt2", db=CONTAMINANT_DBS),
        expand("databases/contaminant/{db}.4.bt2", db=CONTAMINANT_DBS),
        expand("databases/contaminant/{db}.rev.1.bt2", db=CONTAMINANT_DBS),
        expand("databases/contaminant/{db}.rev.2.bt2", db=CONTAMINANT_DBS),

        # expand("databases/functional/{db}{ext}", db=FUNCTIONAL_DBS, ext=['.bck', '.des', '.prj', '.sds', '.ssp', '.suf', '.tis', '-names.txt']),
        # expand("databases/taxonomic/{db}{ext}", db=TAXONOMIC_DBS, ext=['.bck', '.des', '.prj', '.sds', '.ssp', '.suf', '.tis', '-names.txt'])

        # expand("input/{eid}/{sample}.fastq", eid=EID, sample=SAMPLES)

        # samples
        # expand("results/{eid}/joined/{sample}.extendedFrags.fastq", eid=EID, sample=SAMPLES),
        # expand("results/{eid}/joined/{sample}.hist", eid=EID, sample=SAMPLES),
        # expand("results/{eid}/joined/{sample}.notCombined_1.fastq", eid=EID, sample=SAMPLES),
        # expand("results/{eid}/joined/{sample}.notCombined_2.fastq", eid=EID, sample=SAMPLES)


# rule gunzip:
#     input:
#         zipfile = "input/{eid}/{sample}."


rule build_contaminant_references:
    input:
        contaminant_db = "databases/contaminant/{db}.fa"
    output:
        f1 = "databases/contaminant/{db}.1.bt2",
        f2 = "databases/contaminant/{db}.2.bt2",
        f3 = "databases/contaminant/{db}.3.bt2",
        f4 = "databases/contaminant/{db}.4.bt2",
        r1 = "databases/contaminant/{db}.rev.1.bt2",
        r2 = "databases/contaminant/{db}.rev.2.bt2"
    message:
        "Formatting contaminant databases"
    threads: 1
    shell:
        "bowtie2-build {input.contaminant_db} {input.contaminant_db}"


rule build_functional_databases:
    input:
        functional_db = "databases/functional/{db}.fa"
    output:
        f1 = "databases/functional/{db}.bck",
        f2 = "databases/functional/{db}.des",
        f3 = "databases/functional/{db}.prj",
        f4 = "databases/functional/{db}.sds",
        f5 = "databases/functional/{db}.ssp",
        f6 = "databases/functional/{db}.suf",
        f7 = "databases/functional/{db}.tis",
        f8 = "databases/functional/{db}-names.txt"
    message:
        "Formatting functional databases"
    shell:
        "lastdb+ {input.taxonomic_db} {input.taxonomic_db} -p"


rule build_taxonomic_databases:
    input:
        taxonomic_db = "annotation_dbs/taxonomic/{db}.fa"
    output:
        f1 = "databases/taxonomic/{db}.bck",
        f2 = "databases/taxonomic/{db}.des",
        f3 = "databases/taxonomic/{db}.prj",
        f4 = "databases/taxonomic/{db}.sds",
        f5 = "databases/taxonomic/{db}.ssp",
        f6 = "databases/taxonomic/{db}.suf",
        f7 = "databases/taxonomic/{db}.tis",
        f8 = "databases/taxonomic/{db}-names.txt"
    message:
        "Formatting taxonomic databases"
    shell:
        "lastdb+ {input.taxonomic_db} {input.taxonomic_db}"


# rule join_reads:
#     input:
#         r1 = "data/{eid}/{sample}_R1.fastq",
#         r2 = "data/{eid}/{sample}_R2.fastq"
#     output:
#         joined = "results/{eid}/joined/{sample}.extendedFrags.fastq",
#         hist = "results/{eid}/joined/{sample}.hist",
#         failed_r1 = "results/{eid}/joined/{sample}.notCombined_1.fastq",
#         failed_r2 = "results/{eid}/joined/{sample}.notCombined_2.fastq"
#     message: "Joining reads using `flash`"
#     shadow: "shallow"
#     params:
#         min_overlap = config['merging']['minimum_overlap'],
#         max_overlap = config['merging']['maximum_overlap'],
#         max_mismatch_density = config['merging']['maximum_mismatch_density'],
#         phred_offset = config['phred_offset']
#     log: "results/{eid}/logs/{sample}_flash.log"
#     threads: config['merging']['threads']
#     shell: 
#         """flash {input.r1} {input.r2} --min-overlap {params.min_overlap} \
#            --max-overlap {params.max_overlap} --max-mismatch-density {params.max_mismatch_density} \
#            --phred-offset {params.phred_offset} --output-prefix {wildcards.sample} \
#            --output-directory results/{wildcards.eid}/joined/ --threads {threads}"""


# rule filter_contaminants:
#     input:
#         joined = rules.join_reads.output.joined,
#         r1 = rules.join_reads.output.failed_r1,
#         r2 = rules.join_reads.output.failed_r2,
#         prefix = rules.combine_contaminant_references.output
#     output:
#         "results/{eid}/joined/{sample}_joined_filtered.fastq"
#     message:
#         "Aligning all joined and reads that failed to join as single-end against the contamination reference."
#     shadow:
#         "shallow"
#     threads:
#         config['contamination_filtering']['threads']
#     shell:
#         """bowtie2 --threads {threads} --very-sensitive-local -x {input.prefix} -q -U {input.joined},{input.r1},{input.r2} \
#                   | samtools view -@ {threads} -hf4 \
#                   | samtools sort -@ {threads} -T {wildcards.sample} -o -m 8G - \
#                   | bedtools bamtofastq -i stdin -fq {output}
#         """


# rule trim_reads:
#     input:
#         filtered = rule.filter_contaminants.output,
#         # or
#         # filtered = rule.decon_se3.output
#     output:
#         joined = "results/{eid}/trimmed/{sample}_trimmed_filtered_joined.fastq",
#         R1 = "results/{eid}/trimmed/{sample}_trimmed_filtered_R1.fastq",
#         R2 = "results/{eid}/trimmed/{sample}_trimmed_filtered_R2.fastq"
#     params:
#         single_end = config['?']['SE'],
#         phred_value = config['?']['phred33'],
#         min_length = config['?']['length?']
#     message:
#         "Trimming filtered reads using trimmomatic"
#     log:
#         "results/{eid}/trimmed/{sample}.log"
#     shell:
#         """trimmomatic -Xmx32g {params.single_end} -{params.phred_value} {input.filtered} -trimlog {log} \
#             ILLUMINACLIP:adapters/TruSeq2-SE:2:30:10 LEADING:3 TRAILING:3 \
#             SLIDINGWINDOW:4:15 MINLEN:{params.min_length} {output}
#         """


# rule fastqc_R1:
#     input:
#         rule.trim_reads.output.R1
#     output:
#         "results/{eid}/qc/{sample}_R1.fastq"
#     shell:
#         "fastqc {input} -o {output}"


# rule fastqc_R2:
#     input:
#         rule.trim_reads.output.R2
#     output:
#         "results/{eid}/qc/{sample}_R2.fastq"
#     shell:
#         "fastqc {input} -o {output}"


# rule fastqc_joined:
#     input:
#         rule.trim_reads.output.joined
#     output:
#         "results/{eid}/qc/{sample}_joined.fastq"
#     shell:
#         "fastqc {input} -o {output}"

# #######################################################################################

# ######################
# ## Single end reads ##
# ######################

# #For single end reads (Illumina-SE, PacBio, SOiLD)
# rule decon_se:
#     input:
#         db = 'databases/contaminant/{db}'
#         se = rules.reads.output.se
#         prefix = rules.combine_contaminant_references.output
#     output:
#         se = "results/{eid}/decon/{sample}_se.sam",
#     message:
#         "Performing decontamination with bowtie2"
#     threads:
#         config['contamination_filtering']['threads']

#     shell:
#         "bowtie2 --threads {threads} --very-sensitive -x {input.prefix} -q {input.se}"


# rule decon_se2:
#     input:
#         rules.decon.output.se
#     output:
#         se = "results/{eid}/decon/{sample}_se_decon.sam"
#     message:
#         "Filtering with unaligned reads from Decon"
#     shell:
#         """<PICARD> ALIGNMENT_STATUS=Unaligned I={input.se} > {output.se}
#         """


# rule decon_se3:
#     input:
#         rules.decon2.output
#     output:
#         joined = "results/{eid}/decon/{sample}_se_decon.fastq"
#     message:
#         "Converting SAM to FASTQ"
#     shell:
#         """<SAM_to_FASTQ> I={input.se} F={output.se}"""


# rule trim_reads_se:
#     input:
#         filtered = rule.filter_contaminants.output,
#         picard = ?
#     output:
#         SE = "results/{eid}/trimmed/{sample}_trimmed_filtered_se.fastq",
#     params:
#         single_end = config['?']['SE'],
#         phred_value = config['?']['phred33'],
#         min_length = config['?']['length?']
#     message:
#         "Trimming filtered reads using trimmomatic"
#     log:
#         "results/{eid}/trimmed/{sample}.log"
#     shell:
#         """java -Xmx32g -jar trimmomatic-0.33.jar {params.single_end} -{params.phred_value} {input.filtered} -trimlog {log} \
#             ILLUMINACLIP:adapters/TruSeq2-SE:2:30:10 LEADING:3 TRAILING:3 \
#             SLIDINGWINDOW:4:15 MINLEN:{params.min_length} {output}
#         """


# rule fastqc_se:
#     input:
#         rule.trim_reads.output.se
#     output:
#         "results/{eid}/qc/{sample}_se.fastq"
#     shell:
#         "fastqc {input} -o {output}"

# #################################################################################################

# # For metatranscriptomes only!
# rule remove_rRNAs:
#     input:
#         rule.interleave_reads.output,
#         rule.trim_reads.joined.output
#     output:
#         joined = "results/{eid}/trimmed/{sample}_trimmed_filtered_no-rrna_joined.fastq",
#         R1 = "results/{eid}/trimmed/{sample}_trimmed_filtered_no-rrna_R1.fastq",
#         R2 = "results/{eid}/trimmed/{sample}_trimmed_filtered_no-rrna_R2.fastq"
#     message:
#         "Aligning all joined and reads to remove rRNAs for mRNA de novo transcriptome assembly."
#     shadow:
#         "shallow"
#     threads:
#         config['contamination_filtering']['threads']
#     shell:
#         """bowtie2 --threads {threads} --very-sensitive-local -x {input.prefix} -q -U {input.joined},{input.r1},{input.r2} \
#               | samtools view -@ {threads} -hf4 \
#               | samtools sort -@ {threads} -T {wildcards.sample} -o -m 8G - \
#               | bedtools bamtofastq -i stdin -fq {output}
#         """


# rule interleave_reads:
#     input:
#         r1 = rules.join_reads.output.failed_r1,
#         r2 = rules.join_reads.output.failed_r2
#     output:
#         "results/{eid}/trimmed/{sample}_trimmed_filtered_interleaved.fastq"
#     message:
#         "Interleaving non combined R1 and R2 reads"
#     run:
#         IO.interleave_reads(input.r1, input.r2, output)


# rule merge_joined_interleaved:
#     input:
#         il = rules.interleave_reads.output,
#         joined = rules.join_reads.output.joined
#     output:
#         "results/{eid}/assembly/{sample}_merged.fastq"
#     message:
#         "Merging joined and interleaved reads"
#     run:
#         IO.cat_reads(input.il, input.joined, output)


# rule annotate_reads_rRNA:
#     input:
#         rule.interleave_reads.output,
#         trim_reads = rule.trim_reads.joined.output
#     output:
#         rrna = "results/{eid}/annotation/reads/{sample}_{database}"
#     message:
#         "Annotation of rRNAs in quality controlled reads"
#     params:
#         top_hit = config['lastplus']['top_best_hit'],
#         e_value_cutoff = config['lastplus']['e_value_cutoff'],
#         bit_score_cutoff = config['lastplus']['bit_score_cutoff']
#     threads:
#         config['lastplus']['threads']
#     shell:
#         """lastal+ -P {threads} -K {params.top_hit} -E {params.e_value_cutoff} -S {params.bit_score_cutoff} -o {output} \
#             {input.database} {input.trim_reads}"""


# rule taxonomic_placement_rRNAs_LCA:
#     input:
#         rule.annotate_reads_rRNA.output
#     output:
#         rrna_parsed = "results/{eid}/annotation/reads/{sample}_{database}"
#     message:
#         "Parse rRNA reads and place taxonomy using LCA++"
#     params:
#         # TODO
#     threads:
#         config['lastplus']['threads']
#     shell:


# # will want to change this as we add assemblers
# rule assemble:
#     input:
#         # TODO
#     output:
#         # TODO


# #####################
# ##Metagenomes only##
# ####################

# rule megahit:
#     input:
#         rules.filter_contaminants.output
#     output:
#         "results/{eid}/assembly/{sample}.contigs.fa"
#     params:
#         memory = config['assembly']['memory'],
#         min_count = config['assembly']['minimum_count'],
#         k_min = config['assembly']['kmer_min'],
#         k_max = config['assembly']['kmer_max'],
#         k_step = config['assembly']['kmer_step'],
#         merge_level = config['assembly']['merge_level'],
#         prune_level = config['assembly']['prune_level'],
#         low_local_ratio = config['assembly']['low_local_ratio'],
#         min_contig_len = config['assembly']['minimum_contig_length']
#     message:
#         "Assembling using megahit"
#     threads:
#         config['assembly']['threads']
#     shell:
#         """megahit --num-cpu-threads {threads} --memory {params.memory} --read {input} \
#         --k-min {params.k_min} --k-max {params.k_max} --k-step {params.k_step} \
#         --out-dir results/{wildcards.eid}/assembly --out-prefix {wildcards.sample} \
#         --min-contig-len {params.min_contig_len} --min-count {params.min_count} \
#         --merge-level {params.merge_level} --prune-level {params.prune_level} \
#         --low-local-ratio {params.low_local_ratio}"""


# rule metaspades:
#     input:
#         rules.filter_contaminants.output
#     output:
#         "results/results/{eid}/assembly/{sample}"
#     params:
#         memory = config['assembly']['memory']
#     message:
#         "Assembling using metaspades"
#     threads:
#         config['assembly']['threads']
#     shell:
#         """metaspades.py -s1 {input} -o {output} -t {threads} -m {params.memory}"""

# #############################
# ## Metatranscriptomes only ##
# #############################

# rule trinity:
#     input:
#         # need to replace with reference to rule
#         extendedFrags = 'results/{eid}/trimmed/{sample}_trimmed_extendedFrags.fastq',  # after we fix trimming!
#         interleaved = 'results/{eid}/interleaved/{sample}_trimmed_interleaved.fastq'
#     output:
#         "results/{eid}/assembly/{sample}.contigs.fa"
#     params:
#         seqtype = config['assembly']['seqtype'],  # default fastq
#         read_pairing = config['assembly']['single'],  # default single for extendedFrags
#         memory = config['assembly']['max_memory'],
#         run_as_paired = config['assembly']['run_as_paired']
#     message:
#         "Assembling using trinity"
#     threads:
#         config['assembly']['threads']
#     shell:
#         """Trinity --seqType {params.seqtype} --single {input.extendedFrags}, {input.interleaved} \
#         --run_as_paired --max_memory {params.max_memory} --CPU {threads}"""


# rule rnaspades:
#     input:
#         rules.filter_contaminants.output
#     output:
#         "results/results/{eid}/assembly/{sample}"
#     params:
#         memory = config['assembly']['memory']
#     message:
#         "Assembling using metaspades"
#     threads:
#         config['assembly']['threads']
#     shell:
#     """rnaspades.py -s1 {input} -o {output} -t {threads} -m {params.memory}"""


# ########################
# ## Moleculo data only ##
# ########################

# rule truspades:
#     input:
#         input_dir = 'results/{eid}/trimmed/{sample}_trimmed_barcodes1-384_R1.fastq', 'results/{eid}/trimmed/{sample}_trimmed_barcodes1-384_R2.fastq'
#     output:
#         "results/{eid}/assembly/{sample}.contigs.fa"
#     message:
#         "Assembling moleculo data with truspades"
#     threads:
#         config['assembly']['threads']
#     shell:
#         """truspades.py --input-dir input_dir/ -o output_dir/ -t {threads}"""

# #################################################################################################################################3

# rule length_filter:
#     input:
#         rules.assemble.output
#     output:
#         passing = "results/{eid}/assembly/{sample}_length_pass.fa",
#         fail = "results/{eid}/assembly/{sample}_length_fail.fa"
#     params:
#         min_contig_length = config['assembly']['filtered_contig_length']
#     shell:
#         """python scripts/fastx.py length-filter --min-length {params.min_contig_length} \
#         {input} {output.passing} {output.fail}"""


# rule assembly_stats:
#     input:
#         assembled = rules.assembly.output,
#         filtered = rules.length_filter.output
#     output:
#         assembled = "results/{eid}/assembly/{sample}_length_fail_assembly-stats.txt",
#         filtered = "results/{eid}/assembly/{sample}_length_pass_assembly-stats.txt"
#     message:
#         "Obtaining assembly statistics"
#     shell:
#         """perl scripts/CountFasta.pl {input.assembled} > {output.assembled}
#            perl scripts/CountFasta.pl {input.filtered} > {output.filtered}
#         """


# rule merge_assembly_contigs_step1:
#     input:
#         rules.length_filter.output #>1k contigs
#     output:
#         "results/{eid}/assembly/{sample}_cap3-out.fa"
#     message:
#         "Merging contigs with CAP3"
#     shell:
#         """./cap3 {input} -p 95 > {output}.out &"""


# rule length_filter_long_step1:
#     input:
#         rules.assemble.output
#     output:
#         passing = "results/{eid}/assembly/{sample}_length_pass.fa",
#         fail = "results/{eid}/assembly/{sample}_length_fail.fa"
#     params:
#         min_contig_length = config['assembly']['filtered_contig_length']
#     shell:
#         """python scripts/fastx.py length-filter --min-length {params.min_contig_length} \
#         {input} {output.passing} {output.fail}"""


# rule length_filter_long_step2:
#     input:
#         rules.merge_assembly_contigs_step1.output
#     output:
#         passing = "results/{eid}/assembly/{sample}_length_pass.fa",
#         fail = "results/{eid}/assembly/{sample}_length_fail.fa"
#     params:
#         min_contig_length = config['assembly']['filtered_contig_length']
#     shell:
#         """python scripts/fastx.py length-filter --min-length {params.min_contig_length} \
#         {input} {output.passing} {output.fail}"""


# rule merge_long_contigs:
#     input:
#         s1 = rules.length_filter_long_step1.output
#         s2 = rules.length_filter_long_step2.output
#     output:
#         "results/{eid}/assembly/{sample}_merged_contigs.fa"
#     message:
#         "Merging long contigs"
#     run:
#         IO.cat_reads(input.s1, input.s2, output)


# rule merge_assembly_contigs_formatting:
#     input:
#         rules.merge_long_contigs.output
#     output:
#         "results/{eid}/assembly/{sample}_merged_contigs.afg"
#     message:
#         "Formatting all long contigs with minimus2"
#     shell:
#         "toAmos -s {input} -o {input}.afg"


# rule merge_assembly_contigs:
#     input:
#         rules.merge_assembly_contigs_formatting.output
#     output:
#         "results/{eid}/assembly/{sample}_hybrid_assembly.fa"
#     message:
#         "Merging all long contigs with minimus2"
#     params:
#         overlap = config['OVERLAP']
#         maxtrim = contig['MAXTRIM']
#     shell:
#         """minimus2 {input} -D REFCOUNT=0 -D MINID=99.9 -D OVERLAP={params.overlap} \
#         -D MAXTRIM={params.maxtrim} -D WIGGLE=15 -D CONSERR=0.01"""


# rule assembly_stats_hybrid:
#     input:
#         rules.merge_assembly_contigs.output
#     output:
#         hybrid_assembled = "results/{eid}/assembly/{sample}_hybrid_assembly-stats.txt"
#     message:
#         "Obtaining hybrid assembly statistics"
#     shell:
#         """perl scripts/CountFasta.pl {input.assembled} > {output.assembled}"""


# rule map_to_assembly_db_format:
#     input:
#         assembly_db = "assembly/database/{db}"
#         assembly_db_name = os.path.splitext(assembly_db)
#     output:
#         f1 = "{fasta}.1.bt2",
#         f2 = "{fasta}.2.bt2",
#         f3 = "{fasta}.3.bt2",
#         f4 = "{fasta}.4.bt2",
#         r1 = "{fasta}.rev.1.bt2",
#         r2 = "{fasta}.rev.2.bt2"
#     message:
#         "Formatting assembly for mapping"
#     shell:
#         "bowtie2-build {input.assembly_db} {input.assembly_db_name}"


# rule map_to_assembly:
#     input:
#         joined = rules.join_reads.joined
#         interleaved = rules.interleave_reads.output
#         ji = rules.merge_joined_interleaved.output
#         r1 = rules.trim_reads.output.r1
#         r2 = rules.trim_reads.output.r2
#         contigs = rules.megahit.output
#         prefix = os.path.splitext(contigs)
#     output:
#         joined = 'results/{eid}/coverage/{sample}_joined.sam'
#         interleaved = 'results/{eid}/coverage/{sample}_interleaved.sam'
#         ji = 'results/{eid}/coverage/{sample}_joined_interleaved.sam'
#         r1 = 'results/{eid}/coverage/{sample}_r1.sam'
#         r2 = 'results/{eid}/coverage/{sample}_r2.sam'
#     message:
#         "Mapping to assembly"
#     shell:
#         """bowtie2 -p {threads} -x {input.prefix} --very-sensitive-local \
#         -q -U {input.joined}, {input.interleaved}, {input.ji}, {input.r1}, {input.r2}"""


# rule fgsplus_passed:
#     input:
#         rules.length_filter.output.passing
#     output:
#         "results/{eid}/annotation/orfs/{sample}_length_pass.faa",
#     params:
#         sem = config['annotation']['sequencing_error_model'],
#         memory = config['annotation']['memory']
#     message:
#         "Calling protein-coding ORFS with FGS+"
#     threads:
#         config['annotation']['threads']
#     shell:
#         "FGS+ -s {input} -o {output} -w 1 -t {params.sem} -p {threads} -m {params.memory}"


# rule fgsplus_failed:
#     input:
#         rules.length_filter.output.fail
#     output:
#         "results/{eid}/annotation/orfs/{sample}_length_fail.faa"
#     params:
#         sem = config['annotation']['sequencing_error_model'],
#         memory = config['annotation']['memory']
#     message:
#         "Calling protein-coding ORFS with FGS+"
#     threads:
#         config['annotation']['threads']
#     shell:
#         "FGS+ -s {input} -o {output} -w 1 -t {params.sem} -p {threads} -m {params.memory}"


# rule prodigal_orfs_passed:
#     input:
#         rules.length_filter.output.passing
#     output:
#         prot = "results/{eid}/annotation/orfs/{sample}_length_pass.faa",
#         nuc = "results/{eid}/annotation/orfs/{sample}_length_pass.fasta",
#         gff = "results/{eid}/annotation/orfs/{sample}_length_pass.gff"
#     params:
#         g = config['annotation']['translation_table']
#     message:
#         "Calling protein-coding ORFS with prodigal"
#     shell:
#         "prodigal -i {input} -o {output.gff} -f gff -a {output.prot} -d {output.nuc} -g {params.g} -p meta"


# rule prodigal_orfs_failed:
#     input:
#         rules.length_filter.output.fail
#     output:
#         prot = "results/{eid}/annotation/orfs/{sample}_length_fail.faa",
#         nuc = "results/{eid}/annotation/orfs/{sample}_length_fail.fasta",
#         gff = "results/{eid}/annotation/orfs/{sample}_length_fail.gff"
#     params:
#         g = config['annotation']['translation_table']
#     message:
#         "Calling protein-coding ORFS with prodigal"
#     shell:
#         "prodigal -i {input} -o {output.gff} -f gff -a {output.prot} -d {output.nuc} -g {params.g} -p meta"


# rule maxbin_bins:
#     input:
#         reads = rules.filter_contaminants.output
#         contigs = rules.assemble.output
#     output:
#         bins = "results/{eid}/binning/{sample}.fasta",
#         abundance = "results/{eid}/binning/{sample}.abund1",
#         log = "results/{eid}/binning/{sample}.log",
#         marker = "results/{eid}/binning/{sample}.marker",
#         summary = "results/{eid}/binning/{sample}.summary",
#         tooshort = "results/{eid}/binning/{sample}.tooshort",
#         noclass = "results/{eid}/binning/{sample}.noclass"
#     params:
#         min_contig_len = config['binning']['minimum_contig_length'],
#         max_iteration = config['binning']['maximum_iterations'],
#         prob_threshold = config['binning']['probability_threshold'],
#         markerset = config['binning']['marker_set']
#     threads:
#         config['binning']['threads']
#     message:
#         "Binning genomes with MaxBin2.2"
#     shell:
#         """run_MaxBin.pl -contig {input.contigs} -out {output} -reads {input.reads} \
#         -min_contig_length {params.min_contig_len} -max_iteration {params.max_iteration} \
#         -thread {threads} -markerset {params.markerset}"""


# rule lastplus_orfs
#     input:  # how to do wrap rule or ifelse?
#         fgsplus_orfs = rules.fgsplus.output,
#         prodigal_orfs = rules.prodigal.output,
#         database = rules.format_database.output
#     output:
#         annotation = "results/{eid}/annotation/last/{sample}_{database}"
#     params:
#         top_hit = config['lastplus']['top_best_hit'],
#         e_value_cutoff = config['lastplus']['e_value_cutoff'],
#         bit_score_cutoff = config['lastplus']['bit_score_cutoff']
#     threads:
#         config['lastplus']['threads']
#     message:
#         "Annotation of Protein-coding ORFs with Last+"
#     shell:
#         """lastal+ -P {threads} -K {params.top_hit} -E {params.e_value_cutoff} -S {params.bit_score_cutoff} -o {output} \
#         {input.database} {input.fgsplus_orfs}"""

# rule parse_lastplus_lca
#     input:
#         last = rules.lastplus.output
#     output:
#         annotation = "results/{eid}/annotation/last/{sample}_{database}"
#     params:
#         #fill
#     threads
#         #fill
#     message:
#         "Running last+ parsing and taxonomic assignment LCA+"
#     shell:
#         #fill

# rule generate_gtf
#     input:
#         LCA = rules.lcaparselast.output
#     output:
#         gff = "results/{eid}/annotation/quantification/{sample}.gtf"
#     message:
#         "Generate annotated gtf for read quantification"
#     shell:
#         "python src/generate_gff.py {input} {output}"

#  rule run_counting
#     input:
#         gff = rules.generategtf.output
#         sam = rules.maptoassembly.output #aligned only
#     output:
#         verse_counts = "results/{eid}/annotation/quantification/{sample}.counts"
#     message:
#         "Generate counts from reads aligning to contig assembly using VERSE"
#     shell:
#         """verse -a {input.gft} -t {independentassign.default} -g {contig_id} -z 1 -o {output} {input.sam}"""

#  rule get_read_counts #from jeremy zucker scripts
#     input:
#         verse_counts = rules.runcounting.output
#     output:
#         verse_read_counts = "results/{eid}/annotation/quantification/{sample}_read_counts.tsv"
#     message:
#         "Generate read frequencies from VERSE"
#     shell:
#         """python get_read_counts.py --read_count_files {input} --out {output}"""

#  rule get_function_counts #from jeremy zucker scripts
#     input:
#         read_counts = rules.getreadcounts.output
#         annotation_summary_table = rules.lcaparselast.output
#     output:
#         functional_counts_ec = "results/{eid}/annotation/quantification/{sample}_ec.tsv"
#         functional_counts_cog = "results/{eid}/annotation/quantification/{sample}_cog.tsv"
#         functional_counts_ko = "results/{eid}/annotation/quantification/{sample}_ko.tsv"
#         functional_counts_dbcan = "results/{eid}/annotation/quantification/{sample}_dbcan.tsv"
#         functional_counts_metacyc = "results/{eid}/annotation/quantification/{sample}_metacyc.tsv"
#     shell:
#         """python src/get_function_counts.py {input.annotation_summary_table} {input.read_counts} \
#            --function {wildcards.function} --out {output}"""

#  rule get_taxa_counts: #from jeremy zucker scripts
#     input:
#         annotation_summary_table = rules.lcaparselast.output
#         read_counts = rules.getreadcounts.output
#         "Taxonomy/ncbi.map",
#         "Taxonomy/nodes.dmp",
#         "Taxonomy/merged.dmp"
#     output:
#         "{sample}_Counts/{sample}_taxa_{rank}_counts.tsv"
#     shell:
#         "python src/get_rank_counts.py {input} --rank {wildcards.rank} --out {output}"

#  rule funtaxa_counts:
#     input:
#         annotation_summary_table = rules.lcaparselast.output
#         read_counts = rules.getreadcounts.output
#         "Taxonomy/ncbi.map",
#         "Taxonomy/nodes.dmp",
#         "Taxonomy/merged.dmp"
#     output:
#         "{sample}_Counts/{sample}_funtaxa_{function}_{rank}_counts.tsv"
#     shell:
#         """python src/get_fun_rank_counts.py {input} --rank {wildcards.rank} --function {wildcards.function} --out {output}"""
