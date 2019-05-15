import argparse
import os,sys
import pandas as pd
f = open(os.devnull, 'w'); sys.stdout = f # block cufflinks to plot strange code
import plotly.graph_objs as go
from plotly import offline
from cufflinks import iplot
from snakemake.utils import report

PLOTLY_PARAMS = dict(
    include_plotlyjs=False, show_link=False, output_type="div", image_height=700
)

atlas_dir= os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))
sys.path.append(os.path.join(atlas_dir,'scripts'))
from utils.parsers_checkm import read_checkm_output


def main(samples, completeness_files, taxonomy_files, report_out, bin_table):
    sample_data = {}
    div = {}

    df= pd.DataFrame()

    for i,sample in enumerate(samples):
        sample_data= read_checkm_output(taxonomy_table=taxonomy_files[i],
                                        completness_table=completeness_files[i])
        sample_data['Sample']= sample

        df= df.append(sample_data)


    df.to_csv(bin_table,sep='\t')

    div["bin_scatter"] = offline.plot(
        {
            "data": [
                {
                    "x": df.loc[df["Sample"] == sample,"Completeness"],
                    "y": df.loc[df["Sample"] == sample,"Contamination"],
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

.. _ATLAS: https://github.com/metagenome-atlas/atlas

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

    try:
        main(
            samples=snakemake.params.samples,
            taxonomy_files=snakemake.input.taxonomy_files,
            completeness_files=snakemake.input.completeness_files,
            report_out=snakemake.output.report,
            bin_table=snakemake.output.bin_table
        )

    except NameError:
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
