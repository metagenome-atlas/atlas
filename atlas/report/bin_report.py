#! /bin/env python
# handle matplotlib backend
import matplotlib.pylab as plt
plt.switch_backend('agg')
import cufflinks as cf
import os
import pandas as pd
import plotly.graph_objs as go
import numpy as np
import seaborn as sns
import shutil
import zipfile
from plotly import offline
from cufflinks import iplot
from snakemake.utils import report


cf.set_config_file(offline=True, world_readable=True, theme='white')


completeness_files= snakemake.input.completeness_files
taxonomy_files= snakemake.input.taxonomy_files

report_out = snakemake.output.report
table_of_bins = snakemake.output.bin_table

SAMPLES = snakemake.params.samples
folder = os.path.abspath(os.path.dirname(__file__))
stylesheet = os.path.join(folder, 'report.css')


plotly_embed_div = \
    """
        .. raw:: html

            <embed>
                <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
                {}
            </embed>
    """.format
div = {}
report_params = {}
plotly_params = dict(include_plotlyjs=False,
                     show_link=False,
                     output_type='div',
                     image_height=700)





checkm =pd.DataFrame()
for i in range(len(completeness_files)):
    
    comp_= pd.read_table(completeness_files[i],index_col=0)[['Completeness','Contamination','Strain heterogeneity']]
    tax_= pd.read_table(taxonomy_files[i],index_col=0)[['# unique markers (of 43)', '# multi-copy', 'Insertion branch UID',
       'Taxonomy (contained)', 'Taxonomy (sister lineage)', 'GC',
       'Genome size (Mbp)', 'Gene count', 'Coding density']]
    
    checkm=checkm.append(pd.concat((comp_,tax_),axis=1))
    

checkm.to_csv(table_of_bins,sep='\t')    

checkm['Sample']=checkm.index.map(lambda s:s.split('.')[0])


dev= plotly_embed_div(offline.plot(
    {
        'data': [
                    {
                        'x': checkm[checkm['Sample']==sample]['Completeness'],
                        'y': checkm[checkm['Sample']==sample]['Contamination'],
                        'name': sample, 
                        'mode': 'markers',
                        'text': checkm.index[checkm['Sample']==sample],
                        
                        'hoverinfo' : 'text',
                        'showlegend' : False
                    } for sample in checkm.Sample.unique()
                ],
        'layout': {
            'xaxis': {'title': 'Completeness'},
            'yaxis': {'title': "Contamination"}
        }
},**plotly_params))



Report_figure = \
    """
    Recovered bins
    ==============

    {dev}

    See full list at T1_

    """


S= checkm[(checkm.Contamination<=5) &(checkm.Completeness>=90)].copy()


S=S[['Completeness','Contamination', 
     'Taxonomy (contained)','Taxonomy (sister lineage)','GC','Genome size (Mbp)', 'Gene count']]


S['Taxonomy (contained)']=S['Taxonomy (contained)'].apply(lambda s: '; '.join(str(s).split(';')[-2:] ))
S['Taxonomy (sister lineage)']=S['Taxonomy (sister lineage)'].apply(lambda s: '; '.join(str(s).split(';')[-1:] ))


with pd.option_context('display.precision', 3):
    
    table="""
<style>
table {{
    border-collapse: collapse;
    width: 100%;
}}

th, td {{
    text-align: left;
    padding: 8px;
}}

tr:nth-child(even){{background-color: #f2f2f2}}
</style>
<div style="overflow-x:auto;">

{}

</div>

""".format(S.to_html())



table= table.replace('\n','\n'+10*' ')


Report_table=\
    """
    Genomes with > 90% completeness and <5% contamination
    -----------------------------------------------------
    
        .. raw:: html

            <embed>
                {table}
            </embed>
            
    More infmation can be found here:
    
    contigs: {{out_dir}}/{{sample}}/genomic_bins/{{bin}}.fasta
    
    genes:   {{out_dir}}/{{sample}}/genomic_bins/checkm/bins/{{bin}}/genes.faa

            
    """


Combined_report = \
    """
    ####################
    ATLAS binning report
    ####################

    """+Report_figure+Report_table
    
report(Combined_report,
    report_out,T1= table_of_bins)