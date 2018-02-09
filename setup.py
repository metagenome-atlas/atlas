import io
from os.path import dirname, join
from setuptools import setup


long_description = open('README.rst').read()


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
    name='pnnl-atlas',
    version=get_version("atlas/__init__.py"),
    url='https://github.com/pnnl/atlas',
    license='MIT',
    author='Joe Brown',
    author_email='joe.brown@pnnl.gov',
    description='ATLAS - a framework for assembly, annotation, and genomic binning of metagenomic and metatranscriptomic data',
    long_description=long_description,
    packages=['atlas'],
    package_data={'': ['atlas/Snakefile',
                       'atlas/rules/assemble.snakefile',
                       'atlas/rules/annotate.snakefile',
                       'atlas/rules/initialize_checkm.py',
                       'atlas/rules/qc.snakefile',
                       'atlas/envs/optional_genome_binning.yaml',
                       'atlas/envs/required_packages.yaml',
                       'atlas/template_config.yaml',
                       'atlas/report/qc_report.py'
                       ]},
    include_package_data=True,
    # install via conda: click, pandas, pyyaml, snakemake
    install_requires=[],
    entry_points={
          'console_scripts': [
              'atlas = atlas.atlas:cli'
          ]
    },
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
)
