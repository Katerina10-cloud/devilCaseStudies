#!/bin/bash
#SBATCH --mem=400GB
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=THIN
#SBATCH --output=out/devil_human.log
#SBATCH --job-name=devil_human

module load R

# export R_LIBS_USER=~/scratch/r_packages_gpu

Rscript scripts/HumanBlood/cpu_scaling_test_devil.R
