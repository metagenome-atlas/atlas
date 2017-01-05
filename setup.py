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
    url='https://github.com/?????/atlas',
    license='MIT',
    author='Joe Brown',
    author_email='joe.brown@pnnl.gov',
    description='',
    long_description='',
    packages=['atlas'],
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
