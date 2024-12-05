#!/bin/bash
#SBATCH --partition=DGX
#SBATCH --mem=400GB
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=256
#SBATCH --output=out/gpu_macaque.log
#SBATCH --job-name=gpu_macaque
#SBATCH --gpus=8
#module load R/4.3.3

source env.sh
export R_LIBS_USER=~/scratch/r_package_dgx/
Rscript scripts/MacaqueBrain/gpu.R
