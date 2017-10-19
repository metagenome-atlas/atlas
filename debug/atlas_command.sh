#!/bin/bash


source atlas_init.sh
#source activate atlas
#source ../SnakeP/requirements_bbsuite.sh
#atlas assemble config_CR_test.yaml |& tee atlas.log


snakemake -s /home/cegg/kieser/Metagenomics/coding/atlas/atlas/Snakefile \
    -d /home/cegg/kieser/Metagenomics/coding/debug/workidir \
    -p -j 8 \
    --configfile config.yaml \
    --nolock \
    --use-conda \
    --config workflow="complete" $@ |& tee atlas.log
#    --cluster-config /home/cegg/kieser/Metagenomics/coding/debug/slurm.json \
#    --cluster "sbatch -p {cluster.partition} -t {cluster.time} -c {cluster.n} --mem={cluster.mem}g --job-name={cluster.jobname} "\


     #--rerun-incomplete 