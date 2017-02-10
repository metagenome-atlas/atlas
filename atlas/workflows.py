import logging
import os
import sys
from atlas.utils import validate_assembly_config
from subprocess import check_call


def get_snakefile():
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Snakefile")
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf


def assemble(config, jobs, out_dir, dryrun, snakemake_args):
    if not validate_assembly_config(config):
        sys.exit("The configuration file is invalid.")

    cmd = ("snakemake -s {snakefile} -d {out_dir} -p -j {jobs} --configfile '{config}' "
           "--nolock --config workflow=complete {args} --{dryrun}").format(snakefile=get_snakefile(),
                                                                    out_dir=out_dir,
                                                                    jobs=jobs,
                                                                    config=config,
                                                                    dryrun="dryrun" if dryrun else "",
                                                                    args=" ".join(snakemake_args))
    logging.info("Executing: " + cmd)
    check_call(cmd, shell=True)


def download(jobs, out_dir, snakemake_args):
    out_dir = os.path.realpath(out_dir)

    cmd = ("snakemake -s {snakefile} -d {parent_dir} -p -j {jobs} --nolock --rerun-incomplete "
           "--config db_dir='{out_dir}' workflow=download -- {args}").format(snakefile=get_snakefile(),
                                                             parent_dir=os.path.dirname(out_dir),
                                                             jobs=jobs,
                                                             out_dir=out_dir,
                                                             args=" ".join(snakemake_args))

    logging.info("Executing: " + cmd)
    check_call(cmd, shell=True)
