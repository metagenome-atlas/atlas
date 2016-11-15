import argparse
import os
import sys
from pysam import FastxFile


def readfx(fastx):
    """FASTX file reader.

    Args:
        fastx (str): file path to fasta or fastq file; supports gzip compressed files

    Yields:
        tuple: The tuple consists of the read name, the sequence, and quality scores if the file
            is a fastq.

    Raises:
        IOError: If `fastx` file does not exist.

    """
    if not os.path.exists(fastx):
        raise IOError(2, "No such file:", fastx)
    fx = ""
    try:
        fx = FastxFile(fastx)
        for f in fx:
            yield f.name, f.sequence, f.quality
    finally:
        if fx:
            fx.close()


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


def length_sort(fasta, out_file):
    """Sort a fasta file by sequence length.

    Args:
        fasta (str): input fasta file as file path
        out_file (str): output fasta file path

    Returns:
        string: output fasta file path

    """
    records = [(name, seq) for name, seq, qual in readfx(fasta)]
    records.sort(key=lambda x: len(x[1]), reverse=True)
    with open(out_file, 'w') as out_fh:
        for name, seq in records:
            print_fasta_record(name, seq, out_fh)
    return out_file


def length_filter(fasta, pass_file, fail_file, min_len=2000, wrap=80):
    """Filter sequences by length and output a FASTA sorted by length.

    Args:
        fasta (str): input fasta file as file path
        pass_file (str): output file of passing length
        fail_file (str): output file of failing length
        min_len (int): minimum length of passing sequence
        wrap (Optional[int]) : line width of fasta sequence; None is supported for no wrapping

    Returns
        tuple: `pass_file` and `fail_file`

    """
    records = [(name, seq) for name, seq, qual in readfx(fasta)]
    records.sort(key=lambda x: len(x[1]), reverse=True)
    with open(pass_file, 'w') as pass_fh, open(fail_file, 'w') as fail_fh:
        for name, seq in records:
            if len(seq) >= min_len:
                print_fasta_record(name, seq, pass_fh, wrap)
            else:
                print_fasta_record(name, seq, fail_fh, wrap)
    return pass_file, fail_file


def trim_fastx(fastx_file, out_file, left=0, right=0):
    """Trim sequence (and quality) based on left and right options.

    Args:
        fastx_file (str): fasta or fastq sequence file path
        out_file (str): output fastx file path
        left (Optional[int]): number of bases to trim from 5' end
        right (Optional[int]): number of bases to trim from 3' end

    Returns:
        string: out file path

    """
    left = 0 if left < 0 else left
    right = 0 if right < 0 else right

    if left == 0 and right == 0:
        shutil.copy(fastx_file, out_file)
        return out_file

    with open(out_file, 'w') as out_fh:
        for name, seq, qual in readfx(fastx_file):
            seq = seq[left:][:-right] if right > 0 else seq[left:]
            if not seq: continue
            if qual:
                qual = qual[left:][:-right] if right > 0 else qual[left:]
                print('@' + name, seq, '+', qual, sep='\n', file=out_fh)
            else:
                print_fasta_record(name, seq, out_fh)

    return out_file


description = """
length-sort        sort fasta by decreasing length
length-filter      filter fasta by minimum length
trim-fastx         trim 3' or 5' ends of fasta or fastq
"""


def main():

    def pfile_exists(parser, arg):
        if not os.path.exists(arg):
            parser.error("The file %s does not exist." % arg)
        return os.path.abspath(arg)

    fmt = argparse.ArgumentDefaultsHelpFormatter
    p = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="subparser_name", title="commands", description=description)

    length_sort_sp = sp.add_parser('length-sort', formatter_class=fmt, description="TODO")
    length_sort_sp.add_argument("fasta", type=lambda x: pfile_exists(p, x), help="fasta file path")
    length_sort_sp.add_argument("out_file", help="length sorted output fasta file path")
    length_filter_sp = sp.add_parser('length-filter', formatter_class=fmt, description="TODO")
    length_filter_sp.add_argument("fasta", type=lambda x: pfile_exists(p, x), help="fasta file path")
    length_filter_sp.add_argument("pass_file", help="fasta of reads passing threshold")
    length_filter_sp.add_argument("fail_file", help="fasta of reads failing threshold")
    length_filter_sp.add_argument("--min-length", type=int, default=1000, help="minimum length of a passing sequence")
    length_filter_sp.add_argument("--wrap", type=int, default=80, help="sequence wrap length")
    trim_fastx_sp = sp.add_parser('trim-fastx', formatter_class=fmt, description="TODO")
    trim_fastx_sp.add_argument("fastx", type=lambda x: pfile_exists(p, x), help="fasta or fastq file path")
    trim_fastx_sp.add_argument("out_file", help="trimmed fasta or fastq file path")
    trim_fastx_sp.add_argument("-l", "--left", type=int, default=0, help="5' trim length")
    trim_fastx_sp.add_argument("-r", "--right", type=int, default=0, help="3' trim length")
    args = p.parse_args()
    if args.subparser_name == "length-sort":
        length_sort(args.fasta, args.out_file)
    elif args.subparser_name == "length-filter":
        length_filter(args.fasta, args.pass_file, args.fail_file, args.min_length, args.wrap)
    elif args.subparser_name == "trim-fastx":
        trim_fastx(args.fastx_file, args.out_file, args.left, args.right)
    else:
        p.print_help()


if __name__ == "__main__":
    main()
