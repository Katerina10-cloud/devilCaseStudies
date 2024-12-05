#!/bin/bash
#SBATCH --mem=400GB
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=THIN
#SBATCH --output=out/cpu_macaque.log
#SBATCH --job-name=cpu_macaque

module load R
Rscript scripts/MacaqueBrain/cpu.R
