import argparse
import os,sys
f = open(os.devnull, 'w'); sys.stdout = f # block cufflinks to plot strange code
import pandas as pd
import plotly.graph_objs as go
from plotly import offline
from cufflinks import iplot
from snakemake.utils import report

PLOTLY_PARAMS = dict(
    include_plotlyjs=False, show_link=False, output_type="div", image_height=700
)

atlas_dir= os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))

sys.path.append(os.path.join(atlas_dir,'scripts'))


from utils.parsers_bbmap import parse_bbmap_log_file




def parse_map_stats(sample_data, out_tsv):
    stats_df = pd.DataFrame()
    for sample in sample_data.keys():
        df = pd.read_csv(sample_data[sample]["contig_stats"],sep='\t')
        assert df.shape[0] == 1, "Assumed only one row in file {}; found {}".format(
            sample_data[sample]["contig_stats"], df.iloc[0]
        )
        df = df.iloc[0]
        df.name = sample
        genes_df = pd.read_csv(sample_data[sample]["gene_table"], index_col=0,sep='\t')
        df["N_Predicted_Genes"] = genes_df.shape[0]
        used_reads,mapped_reads= parse_bbmap_log_file(sample_data[sample]["mapping_log"])
        df["Assembled_Reads"] = mapped_reads
        df["Percent_Assembled_Reads"] = mapped_reads/used_reads *100

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

.. _ATLAS: https://github.com/metagenome-atlas/atlas

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
    report(report_str, report_out, Table_1=combined_stats, stylesheet=os.path.join(atlas_dir,'report', "report.css"))




if __name__ == "__main__":

    try:
        main(
            samples=snakemake.params.samples,
            contig_stats=snakemake.input.contig_stats,
            gene_tables=snakemake.input.gene_tables,
            mapping_logs=snakemake.input.mapping_logs,
            report_out=snakemake.output.report,
            combined_stats=snakemake.output.combined_contig_stats
        )

    except NameError:

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
