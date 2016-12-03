<p align="center">
    <img src=images/logo.jpg width=300 style="border-radius: 25px; border: 2px solid #000000; align: right" />
</p>

# ATLAS

# Install

# Usage

```
snakemake --jobs 24 \
    --configfile config/atlas_config.yaml \
    --config eid=test-data
```

Preparing to run, place FASTQ files into `data/<eid>`:

```
data/<eid>/<sample-name-1>_R1.fastq
data/<eid>/<sample-name-1>_R2.fastq
data/<eid>/<sample-name-2>_R1.fastq
data/<eid>/<sample-name-2>_R2.fastq
```

## Configuration

# Output
