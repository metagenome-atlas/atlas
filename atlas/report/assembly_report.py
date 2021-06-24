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

    logger.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


# Install exception handler
sys.excepthook = handle_exception


import pandas as pd
import plotly.graph_objs as go
from plotly import offline
from snakemake.utils import report


PLOTLY_PARAMS = dict(
    include_plotlyjs=False, show_link=False, output_type="div", image_height=700
)

atlas_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def main(combined_stats, report_out):

    df = pd.read_csv(combined_stats, sep="\t", index_col=0)
    div = {}
    labels = {
        "Percent_Assembled_Reads": "Percent of Assembled Reads",
        "contig_bp": "Total BP",
        "n_contigs": "Contigs (count)",
        "N_Predicted_Genes": "Predicted Genes (count)",
    }
    for variable in [
        "Percent_Assembled_Reads",
        "contig_bp",
        "n_contigs",
        "N_Predicted_Genes",
    ]:
        y_axis_label = labels[variable]
        div[variable] = offline.plot(
            df[variable].iplot(
                asFigure=True,
                kind="bar",
                xTitle="Samples",
                layout=go.Layout(
                    xaxis=dict(tickangle=45), yaxis=dict(title=y_axis_label)
                ),
            ),
            **PLOTLY_PARAMS,
        )
    div["L50"] = offline.plot(
        df[["L50", "L90"]].iplot(
            asFigure=True,
            kind="bar",
            xTitle="Samples",
            layout=go.Layout(xaxis=dict(tickangle=45), yaxis=(dict(title="Bases"))),
        ),
        **PLOTLY_PARAMS,
    )
    report_str = """

.. raw:: html

    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>


=============================================================
ATLAS_ - Assembly Summary
=============================================================

.. _ATLAS: https://github.com/metagenome-atlas/atlas

.. contents::
    :backlinks: none


Summary
-------

Fragmentation
*************

L50/L90 is a measure of how fractionated assemblies are:
50%/ 90% of the assembly is made up of contigs of Length L50/L90 or longer. Sometimes refered to as N50/N90.


.. raw:: html

    {div[L50]}


Assembly Length
***************

.. raw:: html

    {div[contig_bp]}


Number of Contigs
*****************

.. raw:: html

    {div[n_contigs]}


Number of Predicted Genes
*************************

.. raw:: html

    {div[N_Predicted_Genes]}


Percent of Assembled Reads
**************************

.. raw:: html

    {div[Percent_Assembled_Reads]}


For more information see Table_1_


Downloads
---------

"""
    report(
        report_str,
        report_out,
        Table_1=combined_stats,
        stylesheet=os.path.join(atlas_dir, "report", "report.css"),
    )


if __name__ == "__main__":

    try:

        main(
            combined_stats=snakemake.input.combined_contig_stats,
            report_out=snakemake.output.report,
        )

    except NameError:
        import argparse

        p = argparse.ArgumentParser()
        p.add_argument("--report-out")
        p.add_argument("--combined-stats")
        args = p.parse_args()
        main(
            args.combined_stats,
            args.report_out,
        )
