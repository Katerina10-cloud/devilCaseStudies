
rm(list=ls())
pkgs <- c("ggplot2","dplyr","tidyr","tibble","reshape2", "Seurat", "glmGamPoi", "devil", "nebula")
sapply(pkgs, require, character.only = TRUE)
source("utils.R")

set.seed(SEED)

args = commandArgs(trailingOnly=TRUE)

## Input data
data_path <- args[2]
dataset_name <- args[1]

if (!(file.exists(paste0("results/", dataset_name)))) {
  dir.create(paste0("results/", dataset_name))
}

print(dataset_name)
print(data_path)

seurat_file_name <- paste0('results/', dataset_name, '/seurat.RDS')
if (file.exists(seurat_file_name)) {
  seurat_obj <- readRDS(seurat_file_name)
} else {
  input_data <- read_data(dataset_name, data_path)
  seurat_obj <- prep_seurat_object(input_data, NPC=20, cluster_res = .2)
  saveRDS(seurat_obj, seurat_file_name)
}

time <- dplyr::tibble()
m <- 'devil'
for (m in c("devil", "nebula", "glmGamPoi")) {
  s <- Sys.time()
  de_res_total <- perform_analysis(seurat_obj, method = m)
  e <- Sys.time()
  saveRDS(de_res_total, paste0('results/', dataset_name, '/', m, '.RDS'))
  time <- dplyr::bind_rows(time, dplyr::tibble(method = m, delta_time = e - s))
  print(time)
}

saveRDS(time, paste0('results/', dataset_name, '/time.RDS'))
