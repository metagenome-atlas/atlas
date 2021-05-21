import os, sys

log = open(snakemake.log[0], "w")
sys.stderr = log
sys.stdout = log

import numpy as np
import pandas as pd

pd.options.plotting.backend = "plotly"

from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly import offline


import zipfile
from snakemake.utils import report


PLOTLY_PARAMS = dict(
    include_plotlyjs=False, show_link=False, output_type="div", image_height=700
)


def get_stats_from_zips(zips):
    # def get_read_stats(samples, step):
    quality_pe = pd.DataFrame()
    quality_se = pd.DataFrame()
    for zfile in zips:
        zf = zipfile.ZipFile(zfile)
        # local testing files
        sample = zfile.split("/")[0]  # HACK: sample name is first name of path
        # relative path from snakemake
        # sample = zfile.partition(os.path.sep)[0]

        # single end only
        if "boxplot_quality.txt" in zf.namelist():
            with zf.open("boxplot_quality.txt") as f:
                df = pd.read_csv(f, index_col=0, sep="\t")
                quality_se[sample] = df.mean_1
        else:
            if "se/boxplot_quality.txt" in zf.namelist():
                with zf.open("se/boxplot_quality.txt") as f:
                    df = pd.read_csv(f, index_col=0, sep="\t")
                    quality_se[sample] = df.mean_1

            if "pe/boxplot_quality.txt" in zf.namelist():
                with zf.open("pe/boxplot_quality.txt") as f:
                    df = pd.read_csv(f, index_col=0, sep="\t")
                    df.columns = [df.columns, [sample] * df.shape[1]]
                    quality_pe = pd.concat(
                        (quality_pe, df[["mean_1", "mean_2"]]), axis=1
                    )

    return quality_pe, quality_se


def get_pe_read_quality_plot(df, quality_range, colorscale="Viridis", **kwargs):

    N = len(df["mean_1"].columns)
    c = ["hsl(" + str(h) + ",50%" + ",50%)" for h in np.linspace(0, 360, N + 1)]

    fig = make_subplots(rows=1, cols=2, shared_yaxes=True)

    for i, sample in enumerate(df["mean_1"].columns):
        fig.append_trace(
            dict(
                x=df.index,
                y=df["mean_1"][sample].values,
                type="scatter",
                name=sample,
                legendgroup=sample,
                marker=dict(color=c[i]),
            ),
            1,
            1,
        )

        fig.append_trace(
            dict(
                x=df.index,
                y=df["mean_2"][sample].values,
                type="scatter",
                name=sample,
                legendgroup=sample,
                showlegend=False,
                marker=dict(color=c[i]),
            ),
            1,
            2,
        )

    fig["layout"].update(
        yaxis=dict(range=quality_range, autorange=True, title="Average quality score"),
        xaxis1=dict(title="Position forward read"),
        xaxis2=dict(autorange="reversed", title="Position reverse read"),
    )

    return offline.plot(fig, **kwargs, **PLOTLY_PARAMS)

    #
    # df1 = df[["mean_1", "mean_2"]]
    # return offline.plot(
    #     df1.plot(
    #         subplots=True,
    #         shape=(1, 2),
    #         shared_yaxes=True,
    #         kind="line",
    #         layout = go.Layout(
    #             yaxis = dict(range=quality_range, autorange=True, title="Quality score"),
    #             xaxis1 = dict(title='Position forward read'),
    #             xaxis2 = dict(autorange='reversed',title='Position reverse read')
    #                 )
    #     ),
    #
    # )


# ,


def draw_se_read_quality(df, quality_range, **kwargs):
    return offline.plot(
        df.plot(
            kind="line",
            layout=go.Layout(
                yaxis=dict(
                    range=quality_range, autorange=True, title="Average quality score"
                ),
                xaxis=dict(title="Position read"),
            ),
        ),
        **kwargs,
        **PLOTLY_PARAMS,
    )


def main(report_out, read_counts, zipfiles_QC, min_quality, zipfiles_raw=None):
    div = {}

    # N reads / N bases
    df = pd.read_csv(read_counts, index_col=[0, 1], sep="\t")
    for variable in ["Total_Reads", "Total_Bases"]:

        data = df[variable].unstack()[df.loc[df.index[0][0]].index]

        if "clean" in data.columns:
            data.drop("clean", axis=1, inplace=True)

        div[variable] = offline.plot(
            data.plot(
                kind="bar",
                xTitle="Samples",
                yTitle=variable.replace("_", " "),
                layout=go.Layout(xaxis=dict(tickangle=45)),
            ),
            **PLOTLY_PARAMS,
        )

    Report_numbers = """

Total reads per sample
~~~~~~~~~~~~~~~~~~~~~~

.. raw:: html

    {div[Total_Reads]}

{Legend}

Total bases per sample
~~~~~~~~~~~~~~~~~~~~~~
.. raw:: html

    {div[Total_Bases]}

For details see Table Table1_.
"""
    if data.shape[1] > 1:
        Legend = """
============   ===================================
Step           Output
============   ===================================
raw            the input reads
deduplicated   after (optional) deduplication step
filtered       trimmed, quality filtered
qc             final reads, contaminants removed
============   ===================================
"""
    else:
        Legend = ""

    Report_read_quality_qc = """

Reads quality after QC
~~~~~~~~~~~~~~~~~~~~~~
"""

    Quality_pe, Quality_se = get_stats_from_zips(zipfiles_QC)

    max_quality = 1 + np.nanmax((Quality_pe.max().max(), Quality_se.max().max()))
    if Quality_pe.shape[0] > 0:
        div["quality_qc_pe"] = get_pe_read_quality_plot(
            Quality_pe, [min_quality, max_quality]
        )
        Report_read_quality_qc += """
Paired end
**********
.. raw:: html

    {div[quality_qc_pe]}


"""

    if Quality_se.shape[0] > 0:

        if (Quality_se.shape[0] > 0) & (Quality_se.shape[0] > 0):
            Report_read_quality_qc += """
Single end
+++++++++++

Paired end reads that lost their mate during filtering.

"""

        div["quality_qc_se"] = draw_se_read_quality(
            Quality_se, [min_quality, max_quality]
        )
        Report_read_quality_qc += """

.. raw:: html

    {div[quality_qc_se]}

"""

    if zipfiles_raw is None:
        Report_read_quality_raw = ""
    else:

        Report_read_quality_raw = """

Reads quality before QC
~~~~~~~~~~~~~~~~~~~~~~~

.. raw:: html

    {div[quality_raw]}

"""
        Quality_pe, Quality_se = get_stats_from_zips(zipfiles_raw)
        if Quality_pe.shape[0] > 0:
            div["quality_raw"] = get_pe_read_quality_plot(
                Quality_pe, [min_quality, max_quality]
            )
        elif Quality_se.shape[0] > 0:
            div["quality_raw"] = draw_se_read_quality(
                Quality_se, [min_quality, max_quality]
            )
        else:
            raise IndexError()

    report_str = (
        """

.. raw:: html

    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>


=============================================================
ATLAS_ - QC Summary
=============================================================

.. _ATLAS: https://github.com/metagenome-atlas/atlas

.. contents::
    :backlinks: none


Summary
-------


"""
        + Report_numbers
        + Report_read_quality_qc
        + Report_read_quality_raw
        + """

Downloads
---------

"""
    )

    report(
        report_str,
        report_out,
        Table1=read_counts,
        stylesheet=os.path.join(
            os.path.abspath(os.path.dirname(__file__)), "report.css"
        ),
    )


if __name__ == "__main__":

    try:
        main(
            report_out=snakemake.output.report,
            read_counts=snakemake.input.read_counts,
            zipfiles_raw=snakemake.input.zipfiles_raw
            if hasattr(snakemake.input, "zipfiles_raw")
            else None,
            zipfiles_QC=snakemake.input.zipfiles_QC,
            min_quality=snakemake.params.min_quality,
        )
    except NameError:
        import argparse

        p = argparse.ArgumentParser()
        p.add_argument("--report_out")
        p.add_argument("--read_counts")
        p.add_argument("--zipfiles_raw", nargs="+")
        p.add_argument("--zipfiles_QC", nargs="+")
        p.add_argument("--min_quality")
        args = p.parse_args()

        main(**vars(args))
