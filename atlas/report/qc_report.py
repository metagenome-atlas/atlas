#! /bin/env python
read_counts = snakemake.input.read_counts
report_out = snakemake.output.report



import pandas as pd
from plotly import offline
from cufflinks import iplot
from snakemake.utils import report
import cufflinks as cf
import os, shutil
cf.set_config_file(offline=True, world_readable=True, theme='white')

folder= os.path.abspath(os.path.dirname(__file__))
stylesheet = os.path.join(folder,'report.css')

import plotly.graph_objs as go
layout_x_45 = go.Layout(
    xaxis=dict(
        tickangle=45
    )
)




div={}




df= pd.read_table(read_counts,index_col=[0,1])
for variable in df.columns[:2]:

    data= df[variable].unstack()[df.loc[df.index[0][0]].index]

    div[variable]= offline.plot(data.iplot(asFigure=True,kind='bar',
                                           xTitle='Samples',
                                           yTitle=variable.replace('_',' '),layout=layout_x_45),
                                include_plotlyjs=False,
                                  show_link=False,
                                output_type='div',
                                image_height=700)

Total_Reads=div['Total_Reads']
Total_Bases=div['Total_Bases']
report("""
        ATLAS quality control report
        ===================================

        Total reads per sample
        ----------------------
        .. raw:: html

            <embed>
                <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
                {Total_Reads}
            </embed>

        Legend:

        ============   ===========================================
        Step           Output
        ============   ===========================================
        raw            this are the raw reads given to Atlas
        deduplicated   deduplicated reads are removed (optional)
        filtered       low quality reads are filtered out
        clean          contaminants are removed
        QC             For metagenomic samples 16S reads are added again. This are the final QC reads.
        ============   ===========================================

        Total bases per sample
        ----------------------
        .. raw:: html

            <embed>
                <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
                {Total_Bases}
            </embed>



        for details see Table T1_.
        """, report_out, stylesheet= stylesheet ,T1=read_counts)



# remove ref directory which stores the genome used for _decontamination

if os.path.exists('ref'):
    shutil.rmtree('ref')
