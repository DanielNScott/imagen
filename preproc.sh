#!/bin/sh

###########################
## Slurm Config
###########################
#SBATCH -J imagen-preproc
#SBATCH -t 06:00:00
#SBATCH -n 8
#SBATCH --mem=32000
###########################

# Run the batch job.
Rscript -e 'source("./preproc/import_all.r"); source("./preproc/fmri_routines.r"); import_all("clust")' 1>preproc.out 2>>preproc.err

