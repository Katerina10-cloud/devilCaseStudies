#setwd("~/GitHub/cell_types_analysis")
rm(list=ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble","reshape2", "stringr")
sapply(pkgs, require, character.only = TRUE)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(scMayoMap)
library(Seurat)
source("utils.R")

set.seed(12345)

args = commandArgs(trailingOnly=TRUE)

## Input data
data_path <- args[2]
dataset_name <- args[1]
tissue <- args[3]

dataset_name <- "BaronPancreasData"
tissue <- "pancreas"

if (!(file.exists(paste0("results/", dataset_name)))) {
  dir.create(paste0("results/", dataset_name))
}

if (!(file.exists(paste0("plot/", dataset_name)))) {
  dir.create(paste0("plot/", dataset_name))
}

seurat_obj <- readRDS(paste0('results/', dataset_name, '/seurat.RDS'))
computeGroundTruth(seurat_obj) %>% dplyr::arrange(cluster)
seurat_obj$cell_type <- cell_type_names_to_scMayo_names(seurat_obj$cell_type, tissue)
computeGroundTruth(seurat_obj) %>% dplyr::arrange(cluster)

anno <- scMayoMap::scMayoMapDatabase
anno <- lapply(colnames(anno[2:ncol(anno)]), function(ct) {
  dplyr::tibble(Type=ct, Marker = anno$gene[anno[,ct] == 1])
}) %>% do.call('bind_rows', .)
anno <- anno %>%
  dplyr::filter(grepl(tissue, Type)) %>%
  dplyr::mutate(Type = str_replace_all(Type, paste0(tissue, ":"), "")) %>%
  dplyr::mutate(Type = str_replace_all(Type, paste0(" cell"), ""))

m <- 'devil'
lfc_cut <- 1
pval_cut <- .05
n_markers <- 10

de_res <- de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()
if (sum(grepl("ENSG", de_res$name)) == nrow(de_res)) {
  suppressMessages(de_res$name <- mapIds(org.Hs.eg.db, keys=de_res$name,column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
}

cluster_values <- de_res$cluster %>% unique()
remove_genes <- grepl("^ENS", de_res$name)
de_res <- de_res[!remove_genes, ]

top_genes <- de_res %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(lfc > lfc_cut, adj_pval <= pval_cut) %>%
  dplyr::arrange(-lfc) %>%
  dplyr::slice(1:n_markers)


counts <- as.matrix(seurat_obj@assays$RNA$counts)
percentage_tibble <- lapply(unique(seurat_obj$seurat_clusters), function(c) {
  print(c)

  idx_cluster <- which(seurat_obj$seurat_clusters == c)
  idx_others <- which(!(seurat_obj$seurat_clusters == c))

  if (length(idx_others) > length(idx_cluster)) {
    set.seed(007)
    idx_others <- sample(idx_others, length(idx_cluster), replace = F)
  }

  idxs <- c(idx_cluster, idx_others)

  cluster_pct <- counts[, idx_cluster]
  cluster_pct <- rowSums(cluster_pct > 0) / ncol(cluster_pct)

  non_cluster_pct <- counts[, idx_others]
  non_cluster_pct <- rowSums(non_cluster_pct > 0) / ncol(non_cluster_pct)

  dplyr::tibble(name = names(non_cluster_pct), pct.1 =cluster_pct, pct.2 =non_cluster_pct, cluster = c)
}) %>% do.call('bind_rows', .)

scMayoInput <- lapply(unique(top_genes$cluster), function(c) {
  top_genes %>%
    dplyr::filter(cluster == c) %>%
    dplyr::left_join(percentage_tibble %>% dplyr::filter(cluster == c) %>% dplyr::select(!cluster), by='name') %>%
    dplyr::select(name, pval, adj_pval, lfc, pct.1, pct.2) %>%
    dplyr::rename(p_val_adj = adj_pval, gene=name, avg_log2FC=lfc)
}) %>% do.call('bind_rows', .)

scMayoMap::scMayoMap(scMayoInput, anno, pct.cutoff = 0)

rez <- lapply(unique(top_genes$cluster), function(c) {
  gg <- top_genes$name[top_genes$cluster == c]
  table_pred <- anno %>%
    dplyr::filter(Marker %in% gg) %>%
    dplyr::pull(Type) %>%
    table()
  if (length(table_pred) == 0) {
    dplyr::tibble(pred = NULL, score = 0, model=m, cluster=c)
  } else {
    dplyr::tibble(pred = names(table_pred), score = as.vector(table_pred) / sum(table_pred), model=m, cluster=c)
  }
}) %>% do.call('bind_rows', .) %>%
  dplyr::left_join(computeGroundTruth(seurat_obj), by='cluster') %>%
  na.omit() %>%
  dplyr::mutate(ground_truth = paste0(true_cell_type, " (", cluster, ")"))
#dplyr::mutate(ground_truth = cluster, true_cell_type = cluster)

rez <- rez %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(score == max(score)) %>%
  dplyr::mutate(true_score = ifelse(true_cell_type == pred, 1, 0))

acc = sum(rez$true_score) / length(unique(seurat_obj$seurat_clusters))

comparison_tibble <- dplyr::bind_rows(
  comparison_tibble,
  dplyr::tibble(model=m, acc=acc, n_markers=n_markers, pval_cut=pval_cut)
)











# Prep percentage

