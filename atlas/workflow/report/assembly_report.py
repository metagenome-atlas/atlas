import os, sys
import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logging.captureWarnings(True)


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

#### Begining of scripts

from common_report import *

import os, sys
import pandas as pd
import plotly.express as px


labels = {
    "Percent_Assembled_Reads": "Percent of Assembled Reads",
    "contig_bp": "Total BP",
    "n_contigs": "Contigs (count)",
    "N_Predicted_Genes": "Predicted Genes (count)",
    "N50": "N50-number",
    "L50": "N50-length (bp)",
    "N90": "N90-number",
    "L90": "N90-length (bp)",
}


PLOT_PARAMS = dict(labels=labels)


def make_plots(combined_stats):

    ## Make figures with PLOTLY
    # load and rename data
    df = pd.read_csv(combined_stats, sep="\t", index_col=0)
    df.sort_index(ascending=True, inplace=True)
    df.index.name = "Sample"
    df["Sample"] = df.index

    # create plots store in div
    div = {}

    fig = px.strip(df, y="Percent_Assembled_Reads",hover_name="Sample", **PLOT_PARAMS)
    fig.update_yaxes(range=[0, 100])
    div["Percent_Assembled_Reads"] = fig.to_html(**HTML_PARAMS)

    fig = px.strip(df, y="N_Predicted_Genes", hover_name="Sample",**PLOT_PARAMS)
    div["N_Predicted_Genes"] = fig.to_html(**HTML_PARAMS)

    fig = px.scatter(df, y="L50", x="N50", hover_name="Sample", **PLOT_PARAMS)
    div["N50"] = fig.to_html(**HTML_PARAMS)

    fig = px.scatter(df, y="L90", x="N90", hover_name="Sample", **PLOT_PARAMS)
    div["N90"] = fig.to_html(**HTML_PARAMS)

    fig = px.scatter(
        df, y="contig_bp", x="n_contigs", hover_name="Sample", **PLOT_PARAMS
    )
    div["Total"] = fig.to_html(**HTML_PARAMS)

    return div


# main


div = make_plots(combined_stats=snakemake.input.combined_contig_stats)


make_html(
    div=div,
    report_out=snakemake.output.report,
    html_template_file=os.path.join(reports_dir, "template_assembly_report.html"),
)
