rm(list=ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble","reshape2", "Seurat", "glmGamPoi", "devil", "nebula")
sapply(pkgs, require, character.only = TRUE)

setwd("/orfeo/LTS/CDSLab/LT_storage/kdavydzenka/sc_devil/")
source("devilCaseStudies/multiomics_analysis/utils_multiomics.R")

set.seed(12345)

## Input data
dataset_name <- "MuscleRNA"
data_path <- "data/multiomics/seurat_muscl_rna.RDS"


#data_path <- "data/atac_age.RDS"
#dataset_name <- "MuscleATAC"

if (!(file.exists(paste0("results/", dataset_name)))) {
  dir.create(paste0("results/", dataset_name))
}

input_data <- read_data(dataset_name, data_path)
input_data <- prepare_rna_input(input_data)
#input_data <- prepare_atac_input(input_data)

# RNA analysis #
#time <- dplyr::tibble()
m <- 'nebula'
for (m in c("nebula")) {
 # s <- Sys.time()
  de_res <- perform_analysis_rna(input_data, method = m)
  #e <- Sys.time()
  saveRDS(de_res, paste0('results/', dataset_name, '/', m, '_rna', '.RDS'))
  #time <- dplyr::bind_rows(time, dplyr::tibble(method = m, delta_time = e - s))
  #print(time)
}

# ATAC analysis #
#time <- dplyr::tibble()
#m <- 'devil'
#for (m in c("devil", "nebula", "glmGamPoi")) {
  #s <- Sys.time()
  #de_res <- perform_analysis_atac(input_data, method = m)
  #e <- Sys.time()
  #saveRDS(de_res, paste0('results/', dataset_name, '/', m, '_atac', '.RDS'))
  #time <- dplyr::bind_rows(time, dplyr::tibble(method = m, delta_time = e - s))
  #print(time)
#}

