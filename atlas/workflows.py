import os
import sys
from atlas.utils import validate_assembly_config
from subprocess import check_call


def assemble(config, jobs, out_dir):
    if not validate_assembly_config(config):
        sys.exit("The configuration file is invalid.")
    snakefile = os.path.join(os.path.dirname(os.path.join(os.path.dirname(os.path.abspath(__file__)))), "Snakefile")
    if not os.path.exists(snakefile):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % snakefile)
    try:
        check_call(("snakemake -s {snakefile} -d {out_dir} -p -j {jobs} --configfile {config} "
                    "--nolock --config workflow=complete").format(snakefile=snakefile,
                                                                  out_dir=out_dir,
                                                                  jobs=jobs,
                                                                  config=config), shell=True)
    except:
        # the error will be printed in the snakemake log
        pass
