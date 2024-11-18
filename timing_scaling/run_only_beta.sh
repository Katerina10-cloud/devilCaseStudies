#!/bin/bash
#SBATCH --exclusive
#SBATCH --partition=GPU
#SBATCH --mem=200GB
#BSATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --gpus=2
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --output=out_only_beta.log
#SBATCH --job-name=scaling_de

module load R/4.3.3 cutensor

export R_LIBS_USER=~/scratch/r_packages_gpu

Rscript scaling_only_beta.R
