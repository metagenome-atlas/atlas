import io
from os.path import dirname, join
from setuptools import setup


def get_version(relpath):
  """Read version info from a file without importing it"""
  for line in io.open(join(dirname(__file__), relpath), encoding="cp437"):
    if "__version__" in line:
      if '"' in line:
        # __version__ = "0.9"
        return line.split('"')[1]
      elif "'" in line:
        return line.split("'")[1]


setup(
    name='atlas',
    version=get_version("atlas/__init__.py"),
    url='https://github.com/pnnl/atlas',
    license='MIT',
    author='Joe Brown',
    author_email='joe.brown@pnnl.gov',
    description='',
    long_description='',
    packages=['atlas'],
    package_data={'': ['atlas/Snakefile',
                       'atlas/rules/annotation/diamond.snakefile',
                       'atlas/rules/annotation/diamond.snakefile',
                       'atlas/rules/annotation/munging.snakefile',
                       'atlas/rules/annotation/prodigal.snakefile',
                       'atlas/rules/annotation/verse.snakefile',
                       'atlas/rules/assemblers/megahit.snakefile',
                       'atlas/rules/assemblers/spades.snakefile',
                       'atlas/rules/coassembly/coassemble.snakefile',
                       'atlas/rules/initialization/download.snakefile',
                       'atlas/rules/quality_control/contig_filters.snakefile',
                       'atlas/rules/quality_control/decontamination.snakefile',
                       'atlas/rules/quality_control/error_correction.snakefile',
                       'atlas/rules/quality_control/fastq_filter.snakefile',
                       'atlas/rules/quality_control/fastqc.snakefile',
                       'atlas/rules/quality_control/normalization.snakefile',
                       'atlas/rules/reports/sample.snakefile',
                       ]},
    include_package_data=True,
    install_requires=[
        'click',
        'pandas',
        'pyyaml',
    ],
    entry_points={
          'console_scripts': [
              'atlas = atlas.atlas:cli'
          ]
    },
)
