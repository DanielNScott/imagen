#!/bin/sh

###########################
## Slurm Config
###########################
#SBATCH -J imagen-preproc
#SBATCH -t 6:00:00
#SBATCH -n 1
#SBATCH --nodes 1
#SBATCH --cpus-per-task 16
#SBATCH --mem=32000
###########################

# Run the batch job.
Rscript -e "source('./analysis/import_all.r'); source('./analysis/fmri_routines.r'); import_all('clust',${1},${2})" 1>preproc_${1}_${2}.out 2>>preproc_${1}_${2}.err

