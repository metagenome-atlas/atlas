#!/usr/bin/env python
import os, sys
import logging, traceback


logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logging.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


# Install exception handler
sys.excepthook = handle_exception


# begining of script

import datetime
import shutil
import os


timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%X")


def get_read_stats(fraction, params_in):
    "get read stats by running reformat.sh"

    from snakemake.shell import shell

    subfolder = os.path.join(snakemake.params.folder, fraction)
    tmp_file = os.path.join(subfolder, "read_stats.tmp")
    shell(
        f" mkdir -p {subfolder} 2>> {snakemake.log[0]} "
        " ; "
        f" reformat.sh {params_in} "
        f" bhist={subfolder}/base_hist.txt "
        f" qhist={subfolder}/quality_by_pos.txt "
        f" lhist={subfolder}/readlength.txt "
        f" gchist={subfolder}/gc_hist.txt "
        " gcbins=auto "
        f" bqhist={subfolder}/boxplot_quality.txt "
        f" threads={snakemake.threads} "
        " overwrite=true "
        f" -Xmx{snakemake.resources.java_mem}G "
        f" 2> >(tee -a {snakemake.log[0]} {tmp_file} ) "
    )
    content = open(tmp_file).read()
    pos = content.find("Input:")
    if pos == -1:
        raise Exception("Didn't find read number in file:\n\n" + content)
    else:
        content[pos:].split()[1:4]
        # Input:    123 reads   1234 bases
        n_reads, _, n_bases = content[pos:].split()[1:4]

        os.remove(tmp_file)
    return int(n_reads), int(n_bases)


if len(snakemake.input) >= 2:
    n_reads_pe, n_bases_pe = get_read_stats(
        "pe", "in1={0} in2={1}".format(*snakemake.input)
    )

    n_reads_pe = n_reads_pe / 2

    headers = [
        "Sample",
        "Step",
        "Total_Reads",
        "Total_Bases",
        "Reads_pe",
        "Bases_pe",
        "Reads_se",
        "Bases_se",
        "Timestamp",
    ]

    if os.path.exists(snakemake.params.single_end_file):
        n_reads_se, n_bases_se = get_read_stats(
            "se", "in=" + snakemake.params.single_end_file
        )
    else:
        n_reads_se, n_bases_se = 0, 0

    values = [
        n_reads_pe + n_reads_se,
        n_bases_pe + n_bases_se,
        n_reads_pe,
        n_bases_pe,
        n_reads_se,
        n_bases_se,
    ]
else:
    headers = [
        "Sample",
        "Step",
        "Total_Reads",
        "Total_Bases",
        "Reads",
        "Bases",
        "Timestamp",
    ]
    values = 2 * get_read_stats("", "in=" + snakemake.input[0])

with open(snakemake.output.read_counts, "w") as f:
    f.write("\t".join(headers) + "\n")
    f.write(
        "\t".join(
            [snakemake.wildcards.sample, snakemake.wildcards.step]
            + [str(v) for v in values]
            + [timestamp]
        )
        + "\n"
    )

shutil.make_archive(snakemake.params.folder, "zip", snakemake.params.folder)
shutil.rmtree(snakemake.params.folder)
