rm(list=ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble","reshape2", "Seurat", "glmGamPoi", "devil", "nebula")
sapply(pkgs, require, character.only = TRUE)
source("utils.R")

set.seed(SEED)

args = commandArgs(trailingOnly=TRUE)

## Input data
data_path <- args[2]
dataset_name <- args[1]

# dataset_name <- 'liver'
# data_path <- "datasets/liver.rds"

if (!(file.exists(paste0("results/", dataset_name)))) {
  dir.create(paste0("results/", dataset_name))
}

if (!(file.exists(paste0("plot/", dataset_name)))) {
  dir.create(paste0("plot/", dataset_name))
}

input_data <- read_data(dataset_name, data_path)
seurat_obj <- prep_seurat_object(input_data, NPC=20, cluster_res = .2)

umap_plot_seurat <- Seurat::DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = T,
  repel = T) +
  theme_minimal() +
  theme(legend.position = 'none')

umap_plot_labels <- Seurat::DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "cell_type",
  label = T,
  repel = T) +
  theme_minimal() +
  theme(legend.position = 'none')

saveRDS(seurat_obj, paste0('results/', dataset_name, '/seurat.RDS'))

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
