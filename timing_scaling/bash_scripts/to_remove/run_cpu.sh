#!/bin/bash
#SBATCH --mem=1000GB
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=THIN
#SBATCH --output=out/de_cpu.log
#SBATCH --job-name=cpu_all

module load R

# export R_LIBS_USER=~/scratch/r_packages_gpu

Rscript scripts/cpu_scaling_test_devil.R
Rscript scripts/cpu_scaling_test_glm.R
