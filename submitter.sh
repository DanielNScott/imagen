#!/bin/sh

sbatch ./preproc.sh 1 50
sbatch ./preproc.sh 51 100
sbatch ./preproc.sh 101 150
sbatch ./preproc.sh 151 200
sbatch ./preproc.sh 201 250
sbatch ./preproc.sh 251 300
sbatch ./preproc.sh 301 350
sbatch ./preproc.sh 351 396

