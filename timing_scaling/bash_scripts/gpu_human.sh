#!/bin/bash
#SBATCH --partition=DGX
#SBATCH --mem=400GB
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=256
#SBATCH --output=out/gpu_human.log
#SBATCH --job-name=gpu_human
#SBATCH --gpus=8
#module load R/4.3.3

source env.sh
export R_LIBS_USER=~/scratch/r_package_dgx/
Rscript scripts/HumanBlood/gpu.R
