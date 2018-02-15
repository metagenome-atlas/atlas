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
read_counts = snakemake.input.read_counts
report_out = snakemake.output.report
zip_files_ = snakemake.input.zipfiles
# these are all globals so unsure why only this one received caps
SAMPLES = snakemake.params.samples
folder = os.path.abspath(os.path.dirname(__file__))
stylesheet = os.path.join(folder, 'report.css')
min_quality = snakemake.config["preprocess_minimum_base_quality"]
layout_x_45 = go.Layout(xaxis=dict(tickangle=45))
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

# Plotting functions
# Quality
def get_read_stats(samples, step):
    Quality_pe = pd.DataFrame()
    Quality_se = {}

    for sample in samples:
        zip_file = '{sample}/sequence_quality_control/read_stats/{step}.zip'.format(sample=sample, step=step)
        z = zipfile.ZipFile(zip_file, 'r')

        if 'se/boxplot_quality.txt' in z.namelist():
            with z.open('se/boxplot_quality.txt') as f:
                data = pd.read_table(f, index_col=0)
                Quality_se[sample] = data.mean_1

        Quality_se = pd.DataFrame(Quality_se)
        if 'pe/boxplot_quality.txt' in z.namelist():
            with z.open('pe/boxplot_quality.txt') as f:
                data = pd.read_table(f, index_col=0)
                data.columns = [data.columns, [sample]*data.shape[1]]
                Quality_pe = pd.concat((Quality_pe, data[['mean_1', 'mean_2']]), axis=1)
    return Quality_pe, Quality_se


def draw_pe_read_quality(Quality_pe):
    sns.set_style('whitegrid')
    sns.set_context('talk')
    f, ax = plt.subplots(1, 2, sharey=True, gridspec_kw={'wspace':0.1})
    Quality_pe.mean_1.plot(ax=ax[0], legend=False)
    Quality_pe.mean_2.plot(ax=ax[1], legend=False)
    ax[1].invert_xaxis()

    ax[0].set_ylabel('Quality score')
    ax[0].set_ylim([min_quality, max_quality])

    for i in range(2):
        ax[i].set_xlabel('Position')
        ax[i].set_title('Read {}'.format(i + 1))
    # should pass plotly params as kwargs to this function
    return offline.plot_mpl(f, resize=True, **plotly_params)


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
    div[variable] = plotly_embed_div(offline.plot(data.iplot(asFigure=True, kind='bar',
                                                             xTitle='Samples',
                                                             yTitle=variable.replace('_',' '),
                                                             layout=layout_x_45),
                                                  **plotly_params))
Report_numbers =\
    """
    **********************
    Total reads per sample
    **********************

    {div[Total_Reads]}

    ============   ===========================================
    Step           Output
    ============   ===========================================
    raw            The input reads specified in the configuration
    deduplicated   Counts after (optional) deduplication step
    filtered       Low quality reads have been removed
    clean          Contaminants removed
    QC             Passing reads -- for metagenomic samples, 16S reads are returned to sequence pool
    ============   ===========================================

    Total bases per sample
    ----------------------

    {div[Total_Bases]}

    For details see Table T1_.
    """
report_params['T1'] = read_counts
Report_read_quality_qc=\
    """
    **********************
    Reads quality after QC
    **********************
    """

Quality_pe, Quality_se = get_read_stats(SAMPLES, 'QC')
max_quality = 1 + np.nanmax((Quality_pe.max().max(), Quality_se.max().max()))

if Quality_pe.shape[0] > 0:
    div['quality_qc_pe'] = plotly_embed_div(draw_pe_read_quality(Quality_pe))
    Report_read_quality_qc += \
    """
    Paired end
    ----------

    {div[quality_qc_pe]}

    Single end
    ----------

    Paired end reads that lost their mate during filtering.

    """
div['quality_qc_se'] = plotly_embed_div(draw_se_read_quality(Quality_se))
Report_read_quality_qc += "\n    {div[quality_qc_se]}\n"
Report_read_quality_raw = \
    """
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

### Combined report
Footer = \
    """
    An internet connection is required to vizualize the interactive figures.
    """

Combined_report = \
    """
    ############################
    ATLAS quality control report
    ############################

""" + Report_numbers + Report_read_quality_raw + Report_read_quality_qc + Footer

report(Combined_report, report_out, stylesheet=stylesheet, **report_params)

# remove ref directory which stores the genome used for _decontamination
if os.path.exists('ref'):
    shutil.rmtree('ref')
