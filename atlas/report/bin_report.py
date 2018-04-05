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


def parse_checkm_output(sample_data, out_tsv):
    df = pd.DataFrame()
    for sample in sample_data.keys():
        c_df = pd.read_table(sample_data[sample]["completeness"], index_col=0)[
            ["Completeness", "Contamination", "Strain heterogeneity"]
        ]
        t_df = pd.read_table(sample_data[sample]["taxonomy"], index_col=0)[
            [
                "# unique markers (of 43)",
                "# multi-copy",
                "Insertion branch UID",
                "Taxonomy (contained)",
                "Taxonomy (sister lineage)",
                "GC",
                "Genome size (Mbp)",
                "Gene count",
                "Coding density",
            ]
        ]
        df = df.append(pd.concat([c_df, t_df], axis=1))
    df.to_csv(out_tsv, sep="\t")
    df["Sample"] = df.index.map( lambda s: s.split(".")[0])
    return df


def main(samples, completeness_files, taxonomy_files, report_out, bin_table):
    sample_data = {}
    div = {}
    for sample in samples:
        sample_data[sample] = {}
        for completeness_file in completeness_files:
            # underscore version was for simplified local testing
            # if "%s_" % sample in completeness_file:
            if "%s/" % sample in completeness_file:
                sample_data[sample]["completeness"] = completeness_file
        for taxonomy_file in taxonomy_files:
            # if "%s_" % sample in taxonomy_file:
            if "%s/" % sample in taxonomy_file:
                sample_data[sample]["taxonomy"] = taxonomy_file
    df = parse_checkm_output(sample_data, bin_table)
    div["bin_scatter"] = offline.plot(
        {
            "data": [
                {
                    "x": df[df["Sample"] == sample]["Completeness"],
                    "y": df[df["Sample"] == sample]["Contamination"],
                    "name": sample,
                    "mode": "markers",
                    "text": df.index[df["Sample"] == sample],
                    "hoverinfo": "text",
                    "showlegend": True,
                }
                for sample in df.Sample.unique()
            ],
            "layout": {
                "xaxis": {"title": "Completeness"}, "yaxis": {"title": "Contamination"}
            },
        },
        **PLOTLY_PARAMS,
    )
    # subset the checkm stats dataframe
    df = df[(df.Contamination <= 5) & (df.Completeness >= 90)].copy()
    df.index.name = "Bin ID"
    df.reset_index(inplace=True)
    df = df[
        [
            "Bin ID",
            "Completeness",
            "Contamination",
            "Taxonomy (contained)",
            "Taxonomy (sister lineage)",
            "GC",
            "Genome size (Mbp)",
            "Gene count",
        ]
    ]
    df["Taxonomy (contained)"] = df["Taxonomy (contained)"].apply(
        lambda s: "; ".join(str(s).split(";")[-2:])
    )
    df["Taxonomy (sister lineage)"] = df["Taxonomy (sister lineage)"].apply(
        lambda s: "; ".join(str(s).split(";")[-1:])
    )
    with pd.option_context("display.precision", 3):
        div["table"] = df.to_html(index=False).replace('\n', '\n' + 10 * ' ')
    report_str = """

.. raw:: html

    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>


=============================================================
ATLAS_ - Bin Summary
=============================================================

.. _ATLAS: https://github.com/pnnl/atlas

.. contents::
    :backlinks: none


Summary
-------

Recovered Bins
**************

.. raw:: html

    {div[bin_scatter]}

In some cases, percentages can be above 100% (See: `CheckM Issue 107`_).

.. _CheckM Issue 107: https://github.com/Ecogenomics/CheckM/issues/107


Best Bins
*********

Genomes with >90% completeness and <5% contamination:

.. raw:: html

    {div[table]}


See full list at Table_1_.


Downloads
---------

    """
    report(
        report_str,
        report_out,
        Table_1=bin_table,
        stylesheet=os.path.join(
            os.path.abspath(os.path.dirname(__file__)), "report.css"
        ),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--samples", nargs="+")
    p.add_argument("--completeness", nargs="+")
    p.add_argument("--taxonomy", nargs="+")
    p.add_argument("--report-out")
    p.add_argument("--bin-table")
    args = p.parse_args()
    main(
        args.samples, args.completeness, args.taxonomy, args.report_out, args.bin_table
    )
