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


def make_plots(bin_table):

    div = {}

    div["input_file"] = bin_table

    # Prepare data
    df = pd.read_table(bin_table)

    if snakemake.config["bin_quality_asesser"].lower() == "busco":

        df["Bin Id"] = df["Input_file"].str.replace(".fasta", "", regex=False)

        logging.info("No taxonomic information available, use busco Dataset")

        lineage_name = "Dataset"
        hover_data = [
            "Scores_archaea_odb10",
            "Scores_bacteria_odb10",
            "Scores_eukaryota_odb10",
        ]
        size_name = None

    elif snakemake.config["bin_quality_asesser"].lower() == "checkm":

        df = df.join(
            tax2table(df["Taxonomy (contained)"], remove_prefix=True).fillna("NA")
        )

    
        lineage_name = "phylum"
        size_name = "Genome size (Mbp)"
        hover_data = ["genus"]

    elif snakemake.config["bin_quality_asesser"].lower() == "checkm2":



        df["Bin Id"] = df.index


        lineage_name = "Translation_Table_Used"
        hover_data = ["Completeness_Model_Used","Coding_Density","Contig_N50","GC_Content","Additional_Notes"]
        size_name = "Genome_Size"
    else:
        raise Exception(f"bin_quality_asesser in the config file not understood")

    df.index = df["Bin Id"]

    div[
        "QualityScore"
    ] = "<p>Quality score is calculated as: Completeness - 5 x Contamination.</p>"

    # 2D plot
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

    # By Phylum
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


div = make_plots(bin_table=snakemake.input.bin_table)


make_html(
    div=div,
    report_out=snakemake.output.report,
    html_template_file=os.path.join(reports_dir, "template_bin_report.html"),
    wildcards=snakemake.wildcards,
)
