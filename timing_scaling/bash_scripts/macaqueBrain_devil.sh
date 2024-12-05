#!/bin/bash
#SBATCH --mem=400GB
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=THIN
#SBATCH --output=out/devil_macaque.log
#SBATCH --job-name=devil_macaque

module load R

# export R_LIBS_USER=~/scratch/r_packages_gpu

Rscript scripts/MacaqueBrain/cpu_scaling_test_devil.R
