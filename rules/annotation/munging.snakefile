import json
import os
from itertools import groupby


def read_fasta(fh):
    """Fasta iterator.

    Accepts file handle of .fasta and yields name and sequence.

    Args:
        fh (file): Open file handle of .fasta file

    Yields:
        tuple: name, sequence

    >>> import os
    >>> from itertools import groupby
    >>> f = open("test.fasta", 'w')
    >>> f.write("@seq1\nACTG")
    >>> f.close()
    >>> f = open("test.fastq")
    >>> for name, seq in read_fastq(f):
            assert name == "seq1"
            assert seq == "ACTG"
    >>> f.close()
    >>> os.remove("test.fasta")
    """
    for header, group in groupby(fh, lambda line: line[0] == '>'):
        if header:
            line = next(group)
            name = line[1:].strip()
        else:
            seq = ''.join(line.strip() for line in group)
            yield name, seq


def print_fasta_record(name, seq, out_handle=sys.stdout, wrap=80):
    """Print a fasta record accounting for line wraps in the sequence.

    Args:
        name (str): name or header for fasta entry
        seq (str): sequence
        out_handle (Optional): open file handle in which to write or stdout
        wrap (Optional[int]) : line width of fasta sequence; None is supported for no wrapping

    """
    print('>', name, sep='', file=out_handle)
    if wrap:
        for i in range(0, len(seq), wrap):
            print(seq[i:i + wrap], file=out_handle)
    else:
        print(seq, file=out_handle)


def split_fasta(fasta, chunk_size=250000):
    fasta = os.path.expanduser(fasta)
    root, ext = os.path.splitext(fasta)

    file_idx = 0
    with open(fasta) as f:
        for i, (name, seq) in enumerate(read_fasta(f)):
            if i % chunk_size == 0:
                if i == 0:
                    ofh = open("%s_%d%s" % (root, file_idx, ext), "w")
                    print_fasta_record(name, seq, out_handle=ofh)
                else:
                    ofh.close()
                    file_idx += 1
                    ofh = open("%s_%d%s" % (root, file_idx, ext), "w")
                    print_fasta_record(name, seq, out_handle=ofh)
            else:
                print_fasta_record(name, seq, out_handle=ofh)
    ofh.close()


rule split:
    input:
        "{sample}/annotation/orfs/{sample}.faa"
    output:
        temp(dynamic("{sample}/annotation/orfs/{sample}_{n}.faa"))
    params:
        chunk_size = config["annotation"].get("chunk_size", "250000")
    run:
        split_fasta(input, chunk_size=params.chunk_size)


rule merge_alignments:
    input:
        dynamic("{sample}/annotation/{reference}/{sample}_intermediate_{n}.aln")
    output:
        "{sample}/annotation/{reference}/{sample}_hits.tsv"
    shell:
        "{SHPFXS} cat {input} | sort -k1,1 -k12,12rn > {output}"


rule parse_blast:
    input:
        "{sample}/annotation/{reference}/{sample}_hits.tsv"
    output:
        "{sample}/annotation/{reference}/{sample}_assignments.tsv"
    params:
        # subcommand = lambda wc: "refseq" if "refseq" in wc.reference else "eggnog",
        namemap = lambda wc: config["annotation"]["references"][wc.reference]["namemap"],
        treefile = lambda wc: config["annotation"]["references"][wc.reference].get("tree", ""),
        summary_method = lambda wc: config["annotation"]["references"][wc.reference].get("summary_method", "best"),
        aggregation_method = lambda wc: "--aggregation-method %s" % config["annotation"]["references"][wc.reference].get("aggregation_method", "") if "refseq" in wc.reference else "",
        majority_threshold = lambda wc: "--majority-threshold %f" % config["annotation"]["references"][wc.reference].get("majority_threshold", 0.51) if "refseq" in wc.reference else "",
        min_identity = lambda wc: config["annotation"]["references"][wc.reference].get("min_identity", "50"),
        min_bitscore = lambda wc: config["annotation"]["references"][wc.reference].get("min_bitscore", "0"),
        min_length = lambda wc: config["annotation"]["references"][wc.reference].get("min_length", "60"),
        max_evalue = lambda wc: config["annotation"]["references"][wc.reference].get("max_evalue", "0.000001"),
        max_hits = lambda wc: config["annotation"]["references"][wc.reference].get("max_hits", "10"),
        top_fraction = lambda wc: config["annotation"]["references"][wc.reference].get("top_fraction", "0.50")
    shell:
        """{SHPFXS} python scripts/blast2assignment.py {wildcards.reference} \
               --summary-method {params.summary_method} {params.aggregation_method} \
               {params.majority_threshold} --min-identity {params.min_identity} \
               --min-bitscore {params.min_bitscore} --min-length {params.min_length} \
               --max-evalue {params.max_evalue} --max-hits {params.max_hits} \
               --top-fraction {params.top_fraction} {input} {params.namemap} {params.treefile} \
               {output}"""


rule merge_blast:
    input:
        ["{sample}/annotation/%s/{sample}_assignments.tsv" % i for i in list(config["annotation"]["references"].keys())]
    output:
        "{sample}/annotation/{sample}_merged_assignments.tsv"
    shell:
        "{SHPFXS} python scripts/blast2assignment.py merge-tables {input} {output}"


rule aggregate_counts:
    input:
        merged = "{sample}/annotation/{sample}_merged_assignments.tsv",
        counts = "{sample}/annotation/orfs/{sample}.CDS.txt"
    output:
        ["{sample}/count_tables/{sample}_%s.tsv" % i for i in TABLES]
    params:
        prefix = lambda wc: "{sample}/count_tables/{sample}".format(sample=wc.sample),
        combos = json.dumps(config["summary_counts"])
    shell:
        """{SHPFXS} python scripts/blast2assignment.py counts {params.prefix} {input.merged} \
               {input.counts} '{params.combos}'"""
