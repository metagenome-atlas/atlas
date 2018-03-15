import argparse
import os
import pandas as pd
import plotly.graph_objs as go
from plotly import offline
from cufflinks import iplot
from snakemake.utils import report

PLOTLY_PARAMS = dict(
    include_plotlyjs=False, show_link=False, output_type="div", image_height=700
)


def parse_log_file(log_file, keyword, expect_one_value=True):
    content = open(log_file).read()
    pos = content.find(keyword)
    if pos == -1:
        raise Exception("Didn't find {} in file:\n\n{}".format(keyword, log_file))

    else:
        if expect_one_value:
            return content[pos:].split()[2]

        else:
            return content[pos:].split()[2:]


def parse_map_stats(sample_data, out_tsv):
    stats_df = pd.DataFrame()
    for sample in sample_data.keys():
        df = pd.read_table(sample_data[sample]["contig_stats"])
        assert df.shape[0] == 1, "Assumed only one row in file {}; found {}".format(
            sample_data[sample]["contig_stats"], df.iloc[0]
        )
        df = df.iloc[0]
        df.name = sample
        genes_df = pd.read_table(sample_data[sample]["gene_table"], index_col=0)
        df["N_Predicted_Genes"] = genes_df.shape[0]
        df["Assembled_Reads"] = parse_log_file(
            sample_data[sample]["mapping_log"], "Mapped reads"
        )
        df["Percent_Assembled_Reads"] = parse_log_file(
            sample_data[sample]["mapping_log"], "Percent mapped"
        )
        stats_df = stats_df.append(df)
    stats_df = stats_df.loc[:, ~ stats_df.columns.str.startswith("scaf_")]
    stats_df.columns = stats_df.columns.str.replace("ctg_", "")
    stats_df.to_csv(out_tsv, sep="\t")
    return stats_df


def main(samples, contig_stats, gene_tables, mapping_logs, report_out, combined_stats):
    sample_data = {}
    for sample in samples:
        sample_data[sample] = {}
        for c_stat in contig_stats:
            # underscore version was for simplified local testing
            # if "%s_" % sample in c_stat:
            if "%s/" % sample in c_stat:
                sample_data[sample]["contig_stats"] = c_stat
        for g_table in gene_tables:
            # if "%s_" % sample in g_table:
            if "%s/" % sample in g_table:
                sample_data[sample]["gene_table"] = g_table
        for mapping_log in mapping_logs:
            # if "%s_" % sample in mapping_log:
            if "%s/" % sample in mapping_log:
                sample_data[sample]["mapping_log"] = mapping_log
    df = parse_map_stats(sample_data, combined_stats)
    div = {}
    labels = {
        "Percent_Assembled_Reads": "Percent of Assembled Reads",
        "contig_bp": "Total BP",
        "n_contigs": "Contigs (count)",
        "N_Predicted_Genes": "Predicted Genes (count)",
    }
    for variable in [
        "Percent_Assembled_Reads", "contig_bp", "n_contigs", "N_Predicted_Genes"
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
    div["N50"] = offline.plot(
        df[["N50", "N90"]].iplot(
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

.. _ATLAS: https://github.com/pnnl/atlas

.. contents::
    :backlinks: none


Summary
-------

N50
***

.. raw:: html

    {div[N50]}


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
    report(report_str, report_out, Table_1=combined_stats, stylesheet=os.path.join(os.path.abspath(os.path.dirname(__file__)), "report.css"))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--samples", nargs="+")
    p.add_argument("--contig-stats", nargs="+")
    p.add_argument("--gene-tables", nargs="+")
    p.add_argument("--mapping-logs", nargs="+")
    p.add_argument("--report-out")
    p.add_argument("--combined-stats")
    args = p.parse_args()
    main(
        args.samples,
        args.contig_stats,
        args.gene_tables,
        args.mapping_logs,
        args.report_out,
        args.combined_stats,
    )
