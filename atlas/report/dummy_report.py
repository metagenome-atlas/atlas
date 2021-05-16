import os, sys, stat
log=open(snakemake.log[0],"w")
sys.stderr= log
sys.stdout= log

import pandas as pd
from snakemake.utils import report



def main( input_files,report_out):


input_files= '\n'.join(input_files)


    report_str = """


=============================================================
ATLAS_ - Dummy Report
=============================================================

In this alpha version of atlas the repoprts don't work properly.
I hope they will be available soon.


I protected all input files so they don't get deleted and you will be able to re-create the reports onece the beta version is out.

Input files for this report:
{input_files}



"""
    report(report_str, report_out, stylesheet=os.path.join(atlas_dir,'report', "report.css"))




if __name__ == "__main__":

    input_files =[]
    for file in snakemake.input:

        prrint(f"protect file {file}")
        current_file_status = stat.S_IMODE(os.lstat(file).st_mode)
        os.chmod(file, current_file_status & ~stat.S_IEXEC)
        input_files.append(file)

    main(
        report_out=snakemake.output.report,
        input_files= input_files
    )
