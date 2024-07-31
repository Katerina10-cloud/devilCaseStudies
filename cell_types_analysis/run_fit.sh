#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --mem=200GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=out_fit.log
#SBATCH --job-name=de_analysis

module load R/4.2.3

LC_ALL=C.UTF-8 Rscript 01_fit_data.R BaronPancreasData NULL
LC_ALL=C.UTF-8 Rscript 01_fit_data.R pbmc datasets/pbmc.rds
LC_ALL=C.UTF-8 Rscript 01_fit_data.R liver datasets/liver.rds
LC_ALL=C.UTF-8 Rscript 01_fit_data.R BigBloodData datasets/BigBlood.rds
LC_ALL=C.UTF-8 Rscript 01_fit_data.R BigLiverData datasets/BigLiver.rds
