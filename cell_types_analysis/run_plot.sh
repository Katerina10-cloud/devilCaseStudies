#!/bin/bash
#SBATCH --partition=cpuq
#SBATCH --mem=100GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=out_plot.log
#SBATCH --job-name=de_analysis

module load R/4.2.3

#LC_ALL=C.UTF-8 Rscript 03_supp_figure.R BaronPancreasData pancreas
#LC_ALL=C.UTF-8 Rscript 03_supp_figure.R pbmc blood
#LC_ALL=C.UTF-8 Rscript 03_supp_figure.R liver liver
LC_ALL=C.UTF-8 Rscript 03_supp_figure.R BigBloodData blood
LC_ALL=C.UTF-8 Rscript 03_supp_figure.R BigLiverData liver
