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


logfiles = snakemake.input.mapping_log_files
gene_tables = snakemake.input.gene_tables
contig_stats= snakemake.input.contig_stats


report_out = snakemake.output.report
combined_contig_stats = snakemake.output.combined_contig_stats

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




def parse_log_file(log_file,keyword,expect_one_value=True):

    content = open(log_file).read()

    pos = content.find(keyword)
    if pos == -1:
        raise Exception("Didn't find {} in file:\n\n{}".format(keyword, content))
    else:
        if expect_one_value:
            return content[pos:].split()[2]
        else:
            return content[pos:].split()[2:]




stats =pd.DataFrame()
for i,s in enumerate(SAMPLES):

    d_ = pd.read_table(contig_stats[i])
    assert d_.shape[0]==1, "assumed only one row in file {}. foud {}".format(contig_stats[i],d_.iloc[0])
    d_= d_.iloc[0]
    d_.name= s


    genes= pd.read_table(gene_tables[i],index_col=0)
    d_['n_predicted_genes']= genes.shape[0]

    d_['Assembled_reads']=parse_log_file(logfiles[i],'Mapped reads')
    d_['Percent_assembled_reads']=parse_log_file(logfiles[i],'Percent mapped')


    stats= stats.append(d_)


stats=stats.loc[:,~stats.columns.str.startswith('scaf_')]
stats.columns= stats.columns.str.replace('ctg_','')
stats.to_csv(combined_contig_stats,sep='\t')




div={}

for variable in ['Percent_assembled_reads','contig_bp','n_contigs','n_predicted_genes']:

    div[variable] = plotly_embed_div(offline.plot(stats[variable].iplot(asFigure=True, kind='bar',
                                                             xTitle='Samples',
                                                             layout=go.Layout(xaxis=dict(tickangle=45),yaxis=dict(title=variable.replace('_',' ')))),
                                                  **plotly_params))




div['N50'] = plotly_embed_div(offline.plot(stats[['N50','N90']].iplot(asFigure=True, kind='bar',
                                                         xTitle='Samples',
                                                         layout=go.Layout(xaxis=dict(tickangle=45))),
                                              **plotly_params))



Combined_report = \
    """
    #####################
    ATLAS assembly report
    #####################


    {div[N50]}

    {div[contig_bp]}

    {div[n_contigs]}

    {div[n_predicted_genes]}

    {div[Percent_assembled_reads]}

    More information T1_

    """

report(Combined_report,
    report_out,T1= combined_contig_stats,stylesheet=stylesheet)
