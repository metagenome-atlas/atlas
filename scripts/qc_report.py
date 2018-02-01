from plotly import offline
from cufflinks import iplot
from snakemake.utils import report
import cufflinks as cf
cf.set_config_file(offline=True, world_readable=True, theme='white')


import plotly.graph_objs as go
layout = go.Layout(
    xaxis=dict(
        tickangle=45
    )
)




div={}

read_counts='../read_counts.tsv'

df= pd.read_table(read_counts,index_col=[0,1])
for variable in df.columns[:2]:

    data= df[variable].unstack()[df.loc[df.index[0][0]].index]
    
    div[variable]= offline.plot(data.iplot(asFigure=True,kind='bar',
                                           xTitle='Samples',
                                           yTitle=variable.replace('_',' '),layout=layout), 
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
            
        Total bases per sample
        ----------------------    
        .. raw:: html
        
            <embed>
                <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
                {Total_Bases}
            </embed>
            
            

        for details see Table T1_.
        """, 'report.html', T1=read_counts)
