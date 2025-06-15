#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --mem=36GB
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --output=out/test.log
#SBATCH --job-name=pseudo-replicates

conda activate r_env
Rscript test.R
