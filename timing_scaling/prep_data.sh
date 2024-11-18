#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --mem=500GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=prepping_data.log
#SBATCH --job-name=prep_data

module load R/4.3.3

Rscript make_datasets.R
