#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --mem=200GB
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --output=out.log
#SBATCH --job-name=de_analysis

module load R

LC_ALL=C.UTF-8 Rscript run_models.R "$1"

