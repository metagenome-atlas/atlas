from continuumio/miniconda3



################## METADATA ######################
LABEL base.image="continuumio/miniconda3"
LABEL version="$version"
LABEL software="metagenome-atlas"
LABEL software.version="2"
LABEL about.summary="Three commands to start analysing your metagenome data"
LABEL about.home="https://github.com/metagenome-atlas/atlas"
LABEL about.documentation="https://metagenome-atlas.rtfd.io"
LABEL license="BSD-3"
LABEL about.tags="metagenomics, annotation, snakemake, assembly, genomic-binning, functional-annotation, taxonomic-classifications"

################## MAINTAINER ######################
MAINTAINER Silas Kieser

# Switch back to root for some install

USER root
RUN export LC_ALL=en_US.UTF-8
RUN export LANG=en_US.UTF-8

# setup miniconda
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --set always_yes true

# download atlas
ENV version="2.3.beta2"
RUN wget https://github.com/metagenome-atlas/atlas/archive/${version}.tar.gz
RUN tar -xzf ${version}.tar.gz && mv atlas-${version} atlas
WORKDIR /atlas


#install metagenome atlas
RUN conda env update -n base --file atlasenv.yml
RUN python setup.py install

# short test
RUN atlas --help
RUN atlas --version

ENV databaseDir="/databases"
ENV WORKING_DIR='/.test/Dryrun'

# Dryrun
RUN atlas init --db-dir $databaseDir --threads 3 -w $WORKING_DIR .test/reads/empty
RUN atlas run all -w $WORKING_DIR --dryrun

# download conda packages
RUN atlas run all -w $WORKING_DIR --create-envs-only
RUN atlas run None "logs/checkm_init.txt" -w $WORKING_DIR

# Download databases
# RUN atlas download --db-dir $databaseDir

RUN chmod a+w $databaseDir

# update atlas
RUN rm -rf atlas
RUN pip uninstall -y metagenome-atlas
RUN git clone https://github.com/metagenome-atlas/atlas.git
RUN cd atlas ; pip install .

# short test
RUN atlas --help
RUN atlas --version

# get example data
RUN git clone https://github.com/metagenome-atlas/example_data.git



ENV WORKING_DIR='/AtlasTest'
# Dryrun
RUN atlas init --db-dir $databaseDir --threads 4 --interleaved-fastq -w $WORKING_DIR example_data/reads/test
RUN atlas run all -w $WORKING_DIR --dryrun
RUN atlas run all genomes/tree/checkm.nwk -w $WORKING_DIR --resources mem=50 java_mem=50 --omit-from download_gtdb add_eggNOG_header

RUN atlas run all genomes/tree/checkm.nwk -w $WORKING_DIR --create-envs-only
RUN rm -r example_data
RUN chmod a+w $databaseDir $databaseDir/*

RUN mkdir /WD
WORKDIR /WD

CMD atlas
