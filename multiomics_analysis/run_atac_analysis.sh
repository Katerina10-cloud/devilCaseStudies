#!/bin/bash
#SBATCH --partition=cpuq
#SBATCH --mem=200GB
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=out.log
#SBATCH --job-name=de_analysis

module load R/4.2.3

LC_ALL=C.UTF-8 Rscript ATAC_analysis.R "$1"

