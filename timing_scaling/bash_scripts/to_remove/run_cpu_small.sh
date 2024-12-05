#!/bin/bash
#SBATCH --mem=20GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=EPYC
#SBATCH --output=out/small_de_cpu.log
#SBATCH --job-name=small_de_cpu

module load R
Rscript scripts/cpu_scaling_test_small.R
