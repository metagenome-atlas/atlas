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

#### Begining of scripts


from common_report import *

import pandas as pd
import plotly.express as px


from utils.taxonomy import  tax2table



def make_plots(bin_table):

    div = {}

    div["input_file"] = bin_table

    # Prepare data

    df = pd.read_table(bin_table)
    df = df.join(tax2table(df["Taxonomy (contained)"], remove_prefix=True).fillna("NA"))

    df["Quality Score"] = df.eval("Completeness - 5* Contamination")

    div[
        "QualityScore"
    ] = "<p>Quality score is calculated as: Completeness - 5 x Contamination.</p>"

    # 2D plot
    fig = px.scatter(
        data_frame=df,
        y="Completeness",
        x="Contamination",
        color="phylum",
        size="Genome size (Mbp)",
        hover_data=["genus"],
        hover_name="Bin Id",
    )
    fig.update_yaxes(range=(50, 102))
    fig.update_xaxes(range=(-0.2, 10.1))
    div["2D"] = fig.to_html(**HTML_PARAMS)

    ## By sample
    fig = px.strip(
        data_frame=df,
        y="Quality Score",
        x="Sample",
        color="phylum",
        hover_data=["genus"],
        hover_name="Bin Id",
    )
    fig.update_yaxes(range=(50, 102))
    div["bySample"] = fig.to_html(**HTML_PARAMS)

    # By Phylum
    fig = px.strip(
        data_frame=df,
        y="Quality Score",
        x="phylum",
        hover_data=["genus"],
        hover_name="Bin Id",
    )
    fig.update_yaxes(range=(50, 102))
    div["byPhylum"] = fig.to_html(**HTML_PARAMS)

    return div


## make html report


def make_html(html_template_file, css_file, report_out, div):


# main



div = make_plots(bin_table=snakemake.input.bin_table)

logging.info(div.keys())

make_html(
    div= div,
    report_out=snakemake.output.report,
    html_template_file = os.path.join(reports_dir,"template_bin_report.html" ),
    wildcards= snakemake.wildcards
)
