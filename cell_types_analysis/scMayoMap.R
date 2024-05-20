setwd("~/Desktop/cell_types_analysis")
rm(list=ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble","reshape2")
sapply(pkgs, require, character.only = TRUE)
library(scMayoMap)
library(Seurat)
source("utils.R")

set.seed(12345)

## Input data
DATASET_NAMES
input_data <- read_data('BaronPancreasData')
seurat_obj <- prep_seurat_object(input_data, NPC = 50, cluster_res = .2)

umap_plot_seurat <- Seurat::DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = T,
  repel = T) +
  theme_minimal()

umap_plot_labels <- Seurat::DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "cell_type",
  label = T,
  repel = T) +
  theme_minimal()

umap_plot_seurat
umap_plot_labels

time <- dplyr::tibble()
for (m in c("devil", "nebula", "glm")) {
  s <- Sys.time()
  de_res_total <- perform_analysis(seurat_obj, method = m)
  e <- Sys.time()
}


input_scMayo <- prepScMayoInput(de_res_total, as.matrix(seurat_obj@assays$RNA$counts), seurat_obj$seurat_clusters,
                                n_markers = 50, lfc_cut = 1, pval_cut = 1e-50, distinct_marker = TRUE)

obj <- scMayoMap(data = input_scMayo, tissue = "pancreas")

plt <- scMayoMap.plot(scMayoMap.object = obj, directory = '~/Desktop/', width = 8, height = 6)

de_res_total$is_marker <- lapply(1:nrow(de_res_total), function(i) {
  rr <- de_res_total[i,]
  rr$name
  rr$cluster
  return(nrow(input_scMayo %>% dplyr::filter(cluster == rr$cluster, gene == rr$name)) > 0)
}) %>% unlist()

de_res_total %>%
  dplyr::mutate(adj_pval = if_else(adj_pval <= 1e-50, 1e-50, adj_pval)) %>%
  ggplot(mapping = aes(x=lfc, y=-log10(adj_pval), col=is_marker)) +
  geom_point() +
  facet_wrap(~cluster)

umap1 <- Seurat::DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "seurat_clusters", label = T,
  repel = T) +
  theme_minimal()

umap2 <- Seurat::DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "cell_type",
  label = T,
  repel = T) +
  theme_minimal()

final_plot <- (umap1 | umap2) / plt


#ggsave("final_plot.pdf", dpi=400, width = 16, height = 10, plot = final_plot)

# function to prepare input for scMayoMap
colnames(de_res_total)
de_res <- de_res_total %>% select(name, pval, adj_pval, lfc, cluster)
