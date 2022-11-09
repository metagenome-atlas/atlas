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

# start
import ete3

T = ete3.Tree(snakemake.input.tree, quoted_node_names=True, format=1)

try:

    T.unroot()
    if len(T) > 2:
        T.set_outgroup(T.get_midpoint_outgroup())

except Exception as e:
    logging.error("Failed to root tree, keep unrooted. Reason was:\n\n" + str(e))


T.write(outfile=snakemake.output.tree)
