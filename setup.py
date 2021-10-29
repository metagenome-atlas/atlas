import io
from os.path import dirname, join
from setuptools import setup


# read the contents of your README file
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


# def get_version(relpath):
#     """Read version info from a file without importing it"""
#     for line in io.open(join(dirname(__file__), relpath), encoding="cp437"):
#         if "__version__" in line:
#             if '"' in line:
#                 # __version__ = "0.9"
#                 return line.split('"')[1]
#             elif "'" in line:
#                 return line.split("'")[1]
#





## Get version from git

"""
Gets the current version number.
If in a git repository, it is the current git tag.
Otherwise it is the one contained in the PKG-INFO file.

To use this script, simply import it in your setup.py file
and use the results of get_version() as your package version:

    from version import *

    setup(
        ...
        version=get_version(),
        ...
    )
"""



from os.path import dirname, isdir, join
import os
import re
import subprocess

version_re = re.compile('^Version: (.+)$', re.M)


def get_version():
    d = dirname(__file__)

    if isdir(join(d, '.git')):
        # Get the version using "git describe".
        cmd = 'git describe --tags --match [0-9]*'.split()
        try:
            version = subprocess.check_output(cmd).decode().strip()
        except subprocess.CalledProcessError:
            print('Unable to get version number from git tags')
            exit(1)

        # PEP 386 compatibility
        if '-' in version:
            version = '.post'.join(version.split('-')[:2])

        # Don't declare a version "dirty" merely because a time stamp has
        # changed. If it is dirty, append a ".dev1" suffix to indicate a
        # development revision after the release.
        with open(os.devnull, 'w') as fd_devnull:
            subprocess.call(['git', 'status'],
                            stdout=fd_devnull, stderr=fd_devnull)

        cmd = 'git diff-index --name-only HEAD'.split()
        try:
            dirty = subprocess.check_output(cmd).decode().strip()
        except subprocess.CalledProcessError:
            print('Unable to get git index status')
            exit(1)

        if dirty != '':
            version += '.dev1'

    else:
        # Extract the version from the PKG-INFO file.
        with open(join(d, 'PKG-INFO')) as f:
            version = version_re.search(f.read()).group(1)

    return version









setup(
    name="metagenome-atlas",
    version=get_version(),
    url="https://github.com/metagenome-atlas/atlas",
    license="BSD-3",
    author="Joe Brown, Silas Kieser",
    author_email="brwnjm@gmail.com, silas.kieser@gmail.com",
    description="ATLAS - workflows for assembly, annotation, and genomic binning of metagenomic and metatranscriptomic data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=["atlas"],
    package_data={
        "": [
            "atlas/*",
        ]
    },
    data_files=[(".", ["README.md", "LICENSE.txt"])],
    include_package_data=True,
    install_requires=[],
    # install via conda: click, pandas, pyyaml, snakemake
    entry_points={"console_scripts": ["atlas = atlas.atlas:cli"]},
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
)
