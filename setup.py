import io
from os.path import dirname, join
from setuptools import setup


# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


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
    name='metagenome-atlas',
    version=get_version("atlas/__init__.py"),
    url='https://github.com/metagenome-atlas/atlas',
    license='BSD-3',
    author='Joe Brown, Silas Kieser',
    author_email='brwnjm@gmail.com, silas.kieser@gmail.com',
    description='ATLAS - workflows for assembly, annotation, and genomic binning of metagenomic and metatranscriptomic data.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=['atlas'],
    package_data={'': [
            "atlas/Snakefile",
            "atlas/template_config.yaml",
            "atlas/rules/annotate.snakefile",
            "atlas/rules/assemble.snakefile",
            "atlas/rules/binning.snakefile",
            "atlas/rules/download.snakefile",
            "atlas/rules/gene_annotation.snakefile",
            "atlas/rules/genecatalog.snakefile",
            "atlas/rules/get_fasta_of_bins.py",
            "atlas/rules/initialize_checkm.py",
            "atlas/rules/qc.snakefile",
            "atlas/rules/scg_blank_diamond.rb",
            "atlas/envs/DASTool.yaml",
            "atlas/envs/canopy.yaml",
            "atlas/envs/cd-hit.yaml",
            "atlas/envs/checkm.yaml",
            "atlas/envs/concoct.yaml",
            "atlas/envs/dRep.yaml",
            "atlas/envs/eggNOG.yaml",
            "atlas/envs/maxbin.yaml",
            "atlas/envs/metabat.yaml",
            "atlas/envs/mmseqs.yaml",
            "atlas/envs/optional_genome_binning.yaml",
            "atlas/envs/prokka.yaml",
            "atlas/envs/report.yaml",
            "atlas/envs/required_packages.yaml",
            "atlas/envs/sequence_utils.yaml",
            "atlas/template_config.yaml",
            "atlas/report/qc_report.py",
            "atlas/report/report.css",
            "atlas/report/assembly_report.py",
            "atlas/report/bin_report.py"
                       ]},
    include_package_data=True,
    # install via conda: click, pandas, pyyaml, snakemake
    install_requires=[
        'ruamel.yaml==0.15.35'
    ],
    entry_points={
          'console_scripts': [
              'atlas = atlas.atlas:cli'
          ]
    },
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
)
