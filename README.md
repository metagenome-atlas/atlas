# ATLAS

# Install

All dependencies are installed via [conda](https://www.continuum.io/downloads) using the [bioconda](https://github.com/bioconda/bioconda-recipes) channel.
The workflow and some dependencies require Python 3.5.

To install:

```
conda install -c bioconda \
    bbmap diamond fastqc megahit prodigal samtools snakemake spades verse
```

Or as an isolated environment using our `environment.yml` file:

```
conda env create -f environment.yml
```

And load and unload that environment using `source activate atlas_env`
and `source deactivate atlas_env`, respectively.

In the future we plan to push `atlas` to Bioconda and PyPI, but currently
to install you will need to download or clone `atlas`.

Then within that environment and from within the 'atlas' folder:

```
cd atlas
python setup.py install
```

Following install, `atlas` should be executable:

```
$ atlas -h
Usage: atlas [OPTIONS] COMMAND [ARGS]...

  ATLAS

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  assemble      assembly workflow
  download      download reference files
  make-config   prepopulate a configuration file with samples and defaults
```


# Getting Started

After installing, one needs to download the required databases and create a
sample configuration file.


## Databases

To download the databases and their respective metadata databases:

```
atlas download -o /databases
```

The downloads use approximately 40 GB of disk space.


Configuration File
------------------

To create a simple configuration file (which you can later edit), run:

```
atlas make-config --database-dir /databases config.yaml /my-data
```

In this case, 'my-data' is a directory containing .fastq files similar to:

```
$ tree /my-data
/my-data/
├── Sample-1_R1.fastq.gz
├── Sample-1_R2.fastq.gz
├── Sample-2_R1.fastq.gz
└── Sample-2_R2.fastq.gz
```

Paired-end, single-end (currently with some caveats), and interleaved paired-end
FASTQs are supported.

**Single-end currently works, but parameters to assemblers still need to be
altered such that they are not input as interleaved paired-end.**


## Assembly

After editing your configuration file and adjusting any additional parameters
we run assemblies across our samples using:

```
atlas assemble config.yaml
```

By default, this will write results into our current working directory across
the total number of CPU cores available.

For more complete documentation and advanced configuration options, see our [docs](https://non-broken-link).
