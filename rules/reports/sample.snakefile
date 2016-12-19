rule sample_report:
    input:
        contig_stats = "results/{eid}/{sample}/%s/stats/final_contig_stats.txt" % ASSEMBLER,
        base_comp = "results/{eid}/{sample}/%s/stats/postfilter_base_composition.txt" % ASSEMBLER,
        css = "resources/report.css"
    output:
        html = "results/{eid}/{sample}/{sample}_readme.html"
    shadow:
        "shallow"
    run:
        import pandas as pd

        # contig stats table
        df = pd.read_csv(input.contig_stats, sep="\t")
        contig_stats_csv = "contig_stats.csv"
        df.to_csv(contig_stats_csv,
                  columns=["n_contigs", "contig_bp", "ctg_N50", "ctg_N90", "ctg_max", "gc_avg"],
                  index=False)

        # read base composition across final contigs
        df = pd.read_csv(input.base_comp, sep="\t")
        base_composition_positions = "['%s']" % "', '".join(map(str, df["#Pos"]))
        base_composition_a = "[%s]" % ", ".join(map(str, df["A"]))
        base_composition_c = "[%s]" % ", ".join(map(str, df["C"]))
        base_composition_g = "[%s]" % ", ".join(map(str, df["G"]))
        base_composition_t = "[%s]" % ", ".join(map(str, df["T"]))
        base_composition_n = "[%s]" % ", ".join(map(str, df["N"]))

        report("""

===========================================================================================
Sample Report - Experiment: {wildcards.eid}; Sample: {wildcards.sample}
===========================================================================================

.. raw:: html

    <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>
    <script src="https://code.highcharts.com/highcharts.js"></script>
    <script src="https://code.highcharts.com/modules/exporting.js"></script>
    <script type="text/javascript">

    $(function () {{
        $('#read_composition').highcharts({{
            title: {{text: 'Read Base Composition by Position'}},
            xAxis: {{title: {{text: "Position"}}, categories: {base_composition_positions}}},
            yAxis: {{min: 0, title: {{text: 'Fraction'}}}},
            tooltip: {{}},
            credits: {{enabled: false}},
            legend: {{layout: 'vertical', align: 'right', verticalAlign: 'middle', borderWidth: 0}},
            plotOptions: {{series: {{ marker: {{ enabled: false }} }}, column: {{pointPadding: 0.2, borderWidth: 0}}}},
            series: [{{name: 'A', data: {base_composition_a}}},
                     {{name: 'C', data: {base_composition_c}}},
                     {{name: 'G', data: {base_composition_g}}},
                     {{name: 'T', data: {base_composition_t}}},
                     {{name: 'N', data: {base_composition_n}}}]
            }});
    }});
    </script>

.. contents:: Contents
    :backlinks: none

Read Summary
------------

.. raw:: html

    <div id="read_composition" style="min-width: 310px; height: 500px; margin: 0 auto"></div>


Contig Summary
--------------

.. csv-table::
    :header-rows: 1
    :file: {contig_stats_csv}


               """, output.html, metadata="Author: " + config.get("author", "ATLAS"),
               stylesheet=input.css, contig_stats=input.contig_stats)
