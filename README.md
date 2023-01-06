# Metagenome-Atlas

[![Version](https://anaconda.org/bioconda/metagenome-atlas/badges/version.svg)](https://anaconda.org/bioconda/metagenome-atlas)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/metagenome-atlas.svg?label=Bioconda )](https://anaconda.org/bioconda/metagenome-atlas)
[![Documentation Status](https://readthedocs.org/projects/metagenome-atlas/badge/?version=latest)](https://metagenome-atlas.readthedocs.io/en/latest/?badge=latest)
[![follow on twitter](https://img.shields.io/twitter/follow/SilasKieser.svg?style=social&label=Follow)](https://twitter.com/search?f=tweets&q=%40SilasKieser%20%23metagenomeAtlas&src=typd)


Metagenome-atlas is a easy-to-use metagenomic pipeline based on snakemake. It handles all steps from QC, Assembly, Binning, to Annotation.

![scheme of workflow](resources/images/atlas_list.png?raw=true)

You can start using atlas with three commands:
```
    mamba install -y -c bioconda -c conda-forge metagenome-atlas=2.13.1
    atlas init --db-dir databases path/to/fastq/files
    atlas run all
```

# Webpage

[metagenome-atlas.github.io](https://metagenome-atlas.github.io/)

# Documentation

https://metagenome-atlas.readthedocs.io/

[Tutorial](https://github.com/metagenome-atlas/Tutorial)

# Citation

> ATLAS: a Snakemake workflow for assembly, annotation, and genomic binning of metagenome sequence data.  
> Kieser, S., Brown, J., Zdobnov, E. M., Trajkovski, M. & McCue, L. A.   
> BMC Bioinformatics 21, 257 (2020).  
> doi: [10.1186/s12859-020-03585-4](https://doi.org/10.1186/s12859-020-03585-4)


# Developpment/Extensions

Here are some ideas I work or want to work on when I have time. If you want to contribute or have some ideas let me know via a feature request issue.

- Optimized MAG recovery (e.g. [Spacegraphcats](https://github.com/spacegraphcats/spacegraphcats))
- Integration of viruses/plasmid that live for now as [extensions](https://github.com/metagenome-atlas/virome_atlas)
- Add statistics and visualisations as in [atlas_analyze](https://github.com/metagenome-atlas/atlas_analyze)
- Implementation of most rules as snakemake wrapper
- Cloud execution
- Update to new Snakemake version and use cool reports.
