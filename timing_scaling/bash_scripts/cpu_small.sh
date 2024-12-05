#!/bin/bash
#SBATCH --mem=20GB
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=THIN
#SBATCH --output=out/small_cpu.log
#SBATCH --job-name=cpu_small

module load R

# export R_LIBS_USER=~/scratch/r_packages_gpu

Rscript scripts/BaronPancreas/cpu.R
