import plotly.io as pio
import os, sys

atlas_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

reports_dir = os.path.join(atlas_dir, "report")

sys.path.append(os.path.join(atlas_dir, "scripts"))


pio.templates.default = "simple_white"
HTML_PARAMS = dict(
    include_plotlyjs=False,
    full_html=False,
)


## make html report


def make_html(
    html_template_file,
    report_out,
    div,
    css_file=os.path.join(reports_dir, "report.css"),
    wildcards={},
):
    html_template = open(html_template_file).read()
    css_content = open(css_file).read()

    html_string = html_template.format(div=div, css_content=css_content, **wildcards)

    with open(report_out, "w") as outf:
        outf.write(html_string)
