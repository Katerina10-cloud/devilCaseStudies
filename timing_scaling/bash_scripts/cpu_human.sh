#!/bin/bash
#SBATCH --mem=400GB
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=THIN
#SBATCH --output=out/cpu_human.log
#SBATCH --job-name=cpu_huma

module load R
Rscript scripts/HumanBlood/cpu.R
