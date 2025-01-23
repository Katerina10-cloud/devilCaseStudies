#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --mem=200GB
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --output=out/bca.log
#SBATCH --job-name=bca

module load R
LC_ALL=C.UTF-8 Rscript run_models.R bca
