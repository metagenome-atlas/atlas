from setuptools import setup
import versioneer  # script in directory

__author__ = "Silas Kieser, Joe Brown"
__copyright__ = "Copyright 2021, Silas Kieser"
__email__ = "silas.kieser@gmail.com, brwnjm@gmail.com"
__license__ = "BSD-3"

# read the contents of your README file
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


setup(
    name="metagenome-atlas",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    url="https://github.com/metagenome-atlas/atlas",
    license=__license__,
    author=__author__,
    author_email=__email__,
    zip_safe=False,
    description="ATLAS - workflows for assembly, annotation, and genomic binning of metagenomic and metatranscriptomic data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=["atlas"],
    package_dir={'atlas': 'cli'},
    package_data={
        "": [
            "cli/*",
            "workflow/*",
        ]
    },
    data_files=[(".", ["README.md", "LICENSE.txt"])],
    include_package_data=True,
    install_requires=[],
    # install via conda: click, pandas, pyyaml, snakemake
    entry_points={"console_scripts": ["atlas = cli.atlas:cli"]},
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
)
