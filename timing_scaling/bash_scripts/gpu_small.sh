#!/bin/bash
#SBATCH --partition=DGX
#SBATCH --mem=40GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=256
#SBATCH --output=out/small_gpu.log
#SBATCH --job-name=gpu_small
#SBATCH --gpus=8
#module load R/4.3.3

source env.sh
export R_LIBS_USER=~/scratch/r_package_dgx/
Rscript scripts/BaronPancreas/gpu.R
