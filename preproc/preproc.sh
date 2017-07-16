#!/bin/sh

###########################
## Slurm Config
###########################
#SBATCH -J imagen-preproc
#SBATCH -t 06:00:00
#SBATCH -n 20
#SBATCH --mem=32000
###########################

# Run the batch job.
Rscript -e 'source("import_all.r"); import_all("clust")' 1>preproc.out 2>>preproc.err

