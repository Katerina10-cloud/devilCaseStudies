#!/bin/bash
#SBATCH --mem=400GB
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=THIN
#SBATCH --output=out/glm_macaque.log
#SBATCH --job-name=glm_macaque

module load R

# export R_LIBS_USER=~/scratch/r_packages_gpu

Rscript scripts/MacaqueBrain/cpu_scaling_test_glm.R
