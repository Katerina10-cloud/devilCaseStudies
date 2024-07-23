#!/bin/bash
#SBATCH --partition=cpuq
#SBATCH --mem=32GB
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --output=install.log

module load R/4.3.1

conda env remove -p /home/g.santacatterina/.cache/R/basilisk/1.12.1/0
conda env remove -p /home/g.santacatterina/.cache/R/basilisk/1.12.1/rdevil/0.1.0/pydevil 

LC_ALL=C.UTF-8 Rscript install_rdevil.R
