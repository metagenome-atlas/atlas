import argparse
import os
import pandas as pd
import plotly.graph_objs as go
import zipfile
from plotly import offline, tools
from cufflinks import iplot
from snakemake.utils import report

PLOTLY_PARAMS = dict(
    include_plotlyjs=False, show_link=False, output_type="div", image_height=700
)
read_counts = snakemake.input.read_counts
report_out = snakemake.output.report
zip_files_ = snakemake.input.zipfiles
SAMPLES = snakemake.params.samples
min_quality = snakemake.config["preprocess_minimum_base_quality"]


def get_stats_from_zips(zips):
    # def get_read_stats(samples, step):
    quality_pe = pd.DataFrame()
    quality_se = {}
    for zfile in zips:
        zf = zipfile.ZipFile(zfile)
        # local testing files
        sample = os.path.basename(zfile).partition("_")[0]
        # relative path from snakemake
        # sample = zfile.partition(os.path.sep)[0]
        if "se/boxplot_quality.txt" in zf.namelist():
            with zf.open("se/boxplot_quality.txt") as f:
                df = pd.read_table(f, index_col=0)
                quality_se[sample] = df.mean_1
        if "pe/boxplot_quality.txt" in zf.namelist():
            with zf.open("pe/boxplot_quality.txt") as f:
                df = pd.read_table(f, index_col=0)
                df.columns = [df.columns, [sample] * df.shape[1]]
                quality_pe = pd.concat((quality_pe, df[["mean_1", "mean_2"]]), axis=1)
    return quality_pe, quality_se


def get_pe_read_quality_plot(df, **kwargs):
    # div[variable] = offline.plot(
    #     df[variable].iplot(
    #         asFigure=True,
    #         kind="bar",
    #         xTitle="Samples",
    #         layout=go.Layout(
    #             xaxis=dict(tickangle=45), yaxis=dict(title=y_axis_label)
    #         ),
    #     ),
    #     **PLOTLY_PARAMS,
    # )
    # cf.subplots([df1.figure(kind='heatmap'), df2.figure(kind='heatmap')]).iplot()
    df1 = df[["mean_1", "mean_2"]]
    return offline.plot(
        df1.iplot(
            subplots=True,
            shape=(1, 2),
            shared_yaxes=True,
            asFigure=True,
            kind="bar",
            xTitle="TODO",
        ),
        **kwargs,
    )


def main():
    div = {}


def draw_se_read_quality(Quality_se):
    f, ax = plt.subplots(1, 1)
    Quality_se.plot(legend=False, ax=ax)
    ax.set_ylabel('Quality score')
    ax.set_xlabel('Postition')
    ax.set_ylim([min_quality, max_quality])
    # again, kwargs should be added to this function
    return offline.plot_mpl(f, resize=True, **plotly_params)


# N reads / N bases
df = pd.read_table(read_counts, index_col=[0, 1])
for variable in df.columns[:2]:
    data = df[variable].unstack()[df.loc[df.index[0][0]].index]
    div[variable] = plotly_embed_div(
        offline.plot(
            data.iplot(
                asFigure=True,
                kind='bar',
                xTitle='Samples',
                yTitle=variable.replace('_', ' '),
                layout=layout_x_45,
            ),
            **plotly_params,
        )
    )
Report_numbers = """
    **********************
    Total reads per sample
    **********************

    {div[Total_Reads]}

    ============   ===================================
    Step           Output
    ============   ===================================
    raw            the input reads
    deduplicated   after (optional) deduplication step
    filtered       trimmed, PhiX filtered
    clean          contaminants removed
    qc             passing reads (metag includes SSU)
    ============   ===================================

    Total bases per sample
    ----------------------

    {div[Total_Bases]}

    For details see Table T1_.
    """
report_params['T1'] = read_counts
Report_read_quality_qc = """
    **********************
    Reads quality after QC
    **********************
    """
Quality_pe, Quality_se = get_read_stats(SAMPLES, 'QC')
max_quality = 1 + np.nanmax((Quality_pe.max().max(), Quality_se.max().max()))
if Quality_pe.shape[0] > 0:
    div['quality_qc_pe'] = plotly_embed_div(draw_pe_read_quality(Quality_pe))
    Report_read_quality_qc += """
    Paired end
    ----------

    {div[quality_qc_pe]}

    Single end
    ----------

    Paired end reads that lost their mate during filtering.

    """
div['quality_qc_se'] = plotly_embed_div(draw_se_read_quality(Quality_se))
Report_read_quality_qc += "\n    {div[quality_qc_se]}\n"
Report_read_quality_raw = """
    ***********************
    Reads quality before QC
    ***********************

    {div[quality_raw]}

    """
Quality_pe, Quality_se = get_read_stats(SAMPLES, 'raw')
if Quality_pe.shape[0] > 0:
    div['quality_raw'] = plotly_embed_div(draw_pe_read_quality(Quality_pe))
elif Quality_se.shape[0] > 0:
    div['quality_raw'] = plotly_embed_div(draw_se_read_quality(Quality_se))
else:
    raise IndexError()

    report_str = """

.. raw:: html

    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>


=============================================================
ATLAS_ - QC Summary
=============================================================

.. _ATLAS: https://github.com/pnnl/atlas

.. contents::
    :backlinks: none


Summary
-------



Downloads
---------

"""
    report(
        report_str,
        report_out,
        Table_1=combined_stats,
        stylesheet=os.path.join(
            os.path.abspath(os.path.dirname(__file__)), "report.css"
        ),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--samples", nargs="+")
    ...
    args = p.parse_args()
    main(
        args.samples,
        args.contig_stats,
        args.gene_tables,
        args.mapping_logs,
        args.report_out,
        args.combined_stats,
    )
