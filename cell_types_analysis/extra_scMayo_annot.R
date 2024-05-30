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

dataset_name <- "liver"
tissue <- "liver"

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

umap_plot_seurat <- Seurat::DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = T,
  repel = T)

umap_plot_labels <- Seurat::DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "cell_type",
  label = T,
  repel = T)

umap_plot_seurat | umap_plot_labels

counts <- as.matrix(seurat_obj@assays$RNA$counts)
percentage_tibble <- lapply(unique(seurat_obj$seurat_clusters), function(c) {
  print(c)

  idx_cluster <- which(seurat_obj$seurat_clusters == c)
  idx_others <- which(!(seurat_obj$seurat_clusters == c))

  if (length(idx_others) > length(idx_cluster)) {
    set.seed(SEED)
    idx_others <- sample(idx_others, length(idx_cluster), replace = F)
  }

  idxs <- c(idx_cluster, idx_others)

  cluster_pct <- counts[, idx_cluster]
  cluster_pct <- rowSums(cluster_pct > 0) / ncol(cluster_pct)

  non_cluster_pct <- counts[, idx_others]
  non_cluster_pct <- rowSums(non_cluster_pct > 0) / ncol(non_cluster_pct)

  dplyr::tibble(name = names(non_cluster_pct), pct.1 =cluster_pct, pct.2 =non_cluster_pct, cluster = c)
}) %>% do.call('bind_rows', .)

if (sum(grepl("ENSG", percentage_tibble$name)) == nrow(percentage_tibble)) {
  suppressMessages(percentage_tibble$name <- mapIds(org.Hs.eg.db, keys=percentage_tibble$name,column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
}

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

comparison_tibble <- dplyr::tibble()
lfc_cut <- 1
n_markers <- 50
pval_cut <- 1e-20
m <- 'devil'
for (m in c("devil", "nebula", "glmGamPoi")) {
  de_res <- de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()
  if (sum(grepl("ENSG", de_res$name)) == nrow(de_res)) {
    suppressMessages(de_res$name <- mapIds(org.Hs.eg.db, keys=de_res$name,column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
  }
  #for (n_markers in c(3,5,10,20,50,100,200,300)) {
  for (n_markers in c(5,10,50,100)) {
    print(n_markers)
    for (pval_cut in c(.05, 1e-5, 1e-20, 1e-50)) {
    #for (pval_cut in c(.05, 1e-5, 1e-10, 1e-20, 1e-30, 1e-40, 1e-50)) {
      print(pval_cut)

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

      scMayoInput <- lapply(unique(top_genes$cluster), function(c) {
        top_genes %>%
          dplyr::filter(cluster == c) %>%
          dplyr::left_join(percentage_tibble %>% dplyr::filter(cluster == c) %>% dplyr::select(!cluster), by='name') %>%
          dplyr::select(name, pval, adj_pval, lfc, pct.1, pct.2) %>%
          dplyr::rename(p_val_adj = adj_pval, gene=name, avg_log2FC=lfc)
      }) %>% do.call('bind_rows', .)

      scMayoRes <- scMayoMap::scMayoMap(scMayoInput, scMayoMap::scMayoMapDatabase, pct.cutoff = 0, tissue = tissue)
      rez <- scMayoRes$markers %>%
        dplyr::left_join(computeGroundTruth(seurat_obj), by='cluster') %>%
        dplyr::mutate(pred = str_replace_all(celltype, " cell", "")) %>%
        dplyr::mutate(score = as.numeric(score))

      rez <- rez %>%
        dplyr::group_by(cluster) %>%
        dplyr::filter(score == max(score)) %>%
        dplyr::mutate(true_score = ifelse(true_cell_type == pred, score, 0))

      acc = sum(rez$true_score) / length(unique(seurat_obj$seurat_clusters))

      comparison_tibble <- dplyr::bind_rows(
        comparison_tibble,
        dplyr::tibble(model=m, acc=acc, n_markers=n_markers, pval_cut=pval_cut)
      )
    }
  }
}

comparison_tibble %>%
  #dplyr::filter(pval_cut > 1e-10) %>%
  dplyr::group_by(model, n_markers) %>%
  dplyr::summarise(mean_acc = mean(acc), sd_acc = sd(acc)) %>%
  ggplot(mapping = aes(x=n_markers, y=mean_acc, ymin=mean_acc-sd_acc, ymax=mean_acc+sd_acc, col=model)) +
  geom_pointrange(position=position_dodge(width=0.1)) +
  geom_line(position=position_dodge(width=0.1)) +
  theme_bw() +
  scale_color_manual(values=method_colors) +
  scale_x_continuous(transform = 'log10') +
  labs(x = "N markers", y = "Classification accuracy", col="Algorithm", fill="Algorithm")

comparison_tibble %>%
  #dplyr::filter(pval_cut <= 1e-10) %>%
  # dplyr::group_by(model, n_markers) %>%
  # dplyr::summarise(mean_acc = mean(acc), sd_acc = sd(acc)) %>%
  ggplot(mapping = aes(x=n_markers, y=acc, col=model)) +
  geom_point(position=position_dodge(width=0.1)) +
  geom_smooth() +
  theme_bw() +
  scale_color_manual(values=method_colors) +
  scale_x_continuous(transform = 'log10') +
  labs(x = "N markers", y = "Classification accuracy", col="Algorithm", fill="Algorithm")







# Prep percentage

