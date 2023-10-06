#! /usr/bin/env python


import sys, os
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


import pyfastx


faa_iterator = pyfastx.Fastx(snakemake.input.faa, format="fasta",comment=True)
fna_iterator = pyfastx.Fastx(snakemake.input.fna, format="fasta",comment=True)


with open(snakemake.output.faa, "w") as out_faa, open(
    snakemake.output.fna, "w"
) as out_fna, open(snakemake.output.short, "w") as out_short:
    for name, seq, comment in fna_iterator:
        protein = next(faa_iterator)

        # include gene and corresponding protein if gene passes length threshold
        # or annotation contains prodigal info that it's complete
        if (len(seq) >= snakemake.params.minlength_nt) or ("partial=00" in comment):
            out_fna.write(f">{name} {comment}\n{seq}\n")
            out_faa.write(">{0} {2}\n{1}\n".format(*protein))

        else:
            out_short.write(">{0} {2}\n{1}\n".format(*protein))
