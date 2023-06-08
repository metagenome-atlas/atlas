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

import pandas as pd
import plotly.express as px


from utils.taxonomy import tax2table


def make_plots(bin_info):
    div = {}

    div["input_file"] = bin_info

    # Prepare data
    df = pd.read_table(bin_info)
    df.index = df["Bin Id"]  # need it also as column

    # add species info
    bin2species = pd.read_table(snakemake.input.bins2species, index_col=0)

    logging.info(df.head())

    logging.info(bin2species.head())

    df = df.join(bin2species)

    logging.info(df.describe())

    # Bin Id  Completeness    completeness_general    Contamination   completeness_specific   completeness_model_used translation_table_used  coding_density  contig_n50      average_gene_length      genome_size     gc_content      total_coding_sequences  additional_notes        quality_score   sample  Ambigious_bases Length_contigs  Length_scaffolds N50     N_contigs       N_scaffolds     logN50
    hover_data = [
        "Completeness_Model_Used",
        "Coding_Density",
        "N50",
        "GC_Content",
    ]
    size_name = "Genome_Size"

    lineage_name = "Species"

    div[
        "QualityScore"
    ] = "<p>Quality score is calculated as: Completeness - 5 x Contamination.</p>"

    # 2D plot

    logging.info("make 2d plot")
    fig = px.scatter(
        data_frame=df,
        y="Completeness",
        x="Contamination",
        color=lineage_name,
        size=size_name,
        hover_data=hover_data,
        hover_name="Bin Id",
    )
    fig.update_yaxes(range=(50, 102))
    fig.update_xaxes(range=(-0.2, 10.1))
    div["2D"] = fig.to_html(**HTML_PARAMS)

    ## By sample
    logging.info("plot  by sample")
    fig = px.strip(
        data_frame=df,
        y="Quality_score",
        x="Sample",
        color=lineage_name,
        hover_data=hover_data,
        hover_name="Bin Id",
    )
    fig.update_yaxes(range=(50, 102))
    div["bySample"] = fig.to_html(**HTML_PARAMS)

    # By species
    logging.info("plot by species")
    fig = px.strip(
        data_frame=df,
        y="Quality_score",
        x=lineage_name,
        hover_data=hover_data,
        hover_name="Bin Id",
    )
    fig.update_yaxes(range=(50, 102))
    div["byPhylum"] = fig.to_html(**HTML_PARAMS)

    return div


# main


div = make_plots(bin_info=snakemake.input.bin_info)


make_html(
    div=div,
    report_out=snakemake.output.report,
    html_template_file=os.path.join(reports_dir, "template_bin_report.html"),
    wildcards=snakemake.wildcards,
)
