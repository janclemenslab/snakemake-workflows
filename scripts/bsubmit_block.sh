#!/bin/sh
NAME=$1
SCRIPTDIR="/scratch/clemens10/snakemake-workflows/scripts"
DATADIR="/scratch/clemens10/"
LOGDIR="/scratch/clemens10/$NAME/log"
snakemake --rerun-incomplete --keep-going --nolock --notemp \
    -pr --verbose \
    --use-conda --conda-prefix ~/snake_envs \
    --use-envmodules \
    --cores 1 \
    --jobs 999 \
    --directory . \
    --jobscript "$SCRIPTDIR/jobscript.sh" \
    --cluster "python $SCRIPTDIR/slurm-submit.py {dependencies} $DATADIR $LOGDIR"
