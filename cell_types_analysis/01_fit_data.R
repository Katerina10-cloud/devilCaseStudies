rm(list=ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble","reshape2", "Seurat", "glmGamPoi", "devil", "nebula", 'hrbrthemes')
sapply(pkgs, require, character.only = TRUE)
source("utils.R")

set.seed(12345)

## Input data
dataset_name <- 'BaronPancreasData'
input_data <- read_data(dataset_name)
seurat_obj <- prep_seurat_object(input_data, NPC = 50, cluster_res = .2)
saveRDS(seurat_obj, paste0('results/', dataset_name, '_seurat.RDS'))

time <- dplyr::tibble()
for (m in c("devil", "nebula", "glmGamPoi")) {
  s <- Sys.time()
  de_res_total <- perform_analysis(seurat_obj, method = m)
  e <- Sys.time()
  saveRDS(de_res_total, paste0('results/', dataset_name, '_', m, '.RDS'))
  time <- dplyr::bind_rows(time, dplyr::tibble(method = m, delta_time = e - s))
  print(time)
}

saveRDS(time, paste0('results/', dataset_name, '_time.RDS'))
