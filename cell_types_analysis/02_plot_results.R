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

dataset_name <- "BigLiverData_new"
tissue <- "liver"

if (!(file.exists(paste0("results/", dataset_name)))) {
  dir.create(paste0("results/", dataset_name))
}

if (!(file.exists(paste0("plot/", dataset_name)))) {
  dir.create(paste0("plot/", dataset_name))
}

seurat_obj <- readRDS(paste0('results/', dataset_name, '/seurat.RDS'))
computeGroundTruth(seurat_obj)
seurat_obj$cell_type <- cell_type_names_to_scMayo_names(seurat_obj$cell_type, tissue)
computeGroundTruth(seurat_obj)

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

d_umap <- dplyr::tibble(
  x = seurat_obj@reductions$umap@cell.embeddings[,1],
  y = seurat_obj@reductions$umap@cell.embeddings[,2],
  donor_id = seurat_obj$donor,
  cell_type = seurat_obj$cell_type,
  cluster = seurat_obj$seurat_clusters,
)

saveRDS(d_umap, paste0("results/", dataset_name, "/umap_tibble.rds"))
ggsave(filename = paste0("plot/", dataset_name, "/umap_cluster.pdf"), dpi=400, plot = umap_plot_seurat, width = 8, height = 8)
ggsave(filename = paste0("plot/", dataset_name, "/umap_cell_types.pdf"), dpi=400, plot = umap_plot_labels, width = 8, height = 8)

m <- 'devil'
lfc_cut <- 1
pval_cut <- .05
n_markers <- 10

anno <- scMayoMap::scMayoMapDatabase
anno <- lapply(colnames(anno[2:ncol(anno)]), function(ct) {
  dplyr::tibble(Type=ct, Marker = anno$gene[anno[,ct] == 1])
}) %>% do.call('bind_rows', .)
anno <- anno %>%
  dplyr::filter(grepl(tissue, Type)) %>%
  dplyr::mutate(Type = str_replace_all(Type, paste0(tissue, ":"), "")) %>%
  dplyr::mutate(Type = str_replace_all(Type, paste0(" cell"), ""))

anno$Type %>% unique()
seurat_obj$cell_type %>% unique()

m <- 'devil'
for (m in c("devil", "nebula", "glmGamPoi")) {
  de_res <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()

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
    # dplyr::filter(abs(lfc) > lfc_cut, adj_pval <= pval_cut) %>%
    # dplyr::arrange(-abs(lfc)) %>%
    dplyr::slice(1:n_markers) %>%
    dplyr::mutate(gene_id = paste(name, cluster, sep=':'))

  c <- unique(top_genes$cluster)[4]
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
    dplyr::mutate(ground_truth = paste0(true_cell_type, " (", cluster, ")")) %>%
    dplyr::filter(score > .1)
    #dplyr::mutate(ground_truth = cluster, true_cell_type = cluster)

  rez$true_score = lapply(1:nrow(rez), function(i) {
    r <- rez[i,]
    if (r$true_cell_type == r$pred) {
      return(r$score)
    } else {
      return(0)
    }
  }) %>% unlist()

  rez_max <- rez %>%
    dplyr::group_by(cluster) %>%
    dplyr::filter(score == max(score)) %>%
    dplyr::mutate(is_correct = true_cell_type == pred)

  res_heatmap <- ggplot() +
    geom_point(rez %>% dplyr::filter(score > .1), mapping = aes(x = ground_truth, y=pred, col=score, size=score), shape=20) +
    geom_point(rez_max %>% dplyr::filter(score > .1, is_correct), mapping = aes(x = ground_truth, y=pred, size=score), shape=21, col="blue") +
    geom_point(rez_max %>% dplyr::filter(score > .1, !is_correct), mapping = aes(x = ground_truth, y=pred, size=score), shape=21, col="red") +
    theme_bw() +
    scale_color_gradient(low = "white", high = "#2ca25f") +
    #scale_color_continuous(type = "viridis") +
    labs(x = "Ground truth", y = "Prediction") +
    theme(axis.text.x = element_text(angle = 45, vjust = .5))

  volc <- de_res %>%
    dplyr::mutate(name_id = paste(name, cluster, sep = ":")) %>%
    dplyr::mutate(is_marker = name_id %in% top_genes$gene_id) %>%
    dplyr::mutate(adj_pval = ifelse(adj_pval <= 6.181186e-16, 6.181186e-16, adj_pval)) %>%
    ggplot(mapping = aes(x=lfc, y=-log10(adj_pval), col=is_marker)) +
    geom_point() +
    scale_color_manual(values = c(alpha("black", .1), alpha("indianred", 1))) +
    xlim(c(-5, NA)) +
    facet_wrap(~paste0("C ",cluster), scales = "free_x", ncol = length(unique(seurat_obj$seurat_clusters))) +
    #geom_hline(yintercept = -log10(.05), linetype = "dashed") +
    theme_minimal() +
    labs(x = "Fold Change (log2)", y = '-Log10 P', col="Marker") #+
  #theme(legend.position = 'bottom')

  # c <- 1
  # gt_c <- computeGroundTruth(seurat_obj) %>% dplyr::filter(cluster == c) %>% pull(true_cell_type)
  # gg_c <- anno %>% dplyr::filter(Type == gt_c) %>% pull(Marker)
  # de_res %>%
  #   dplyr::filter(cluster == c) %>%
  #   dplyr::mutate(is_marker = name %in% gg_c) %>%
  #   ggplot(mapping = aes(x=lfc, y=-log10(adj_pval), col=is_marker)) +
  #   geom_point() +
  #   scale_color_manual(values = c(alpha("black", .1), alpha("indianred", 1)))

  #obj <- scMayoMap(data = input_scMayo, tissue = tissue, pct.cutoff = 0)
  #saveRDS(obj, paste0("results/", dataset_name, "/", m, "_scMayo.rds"))

  #mayoMatrix <- plot_scMayoOutput(obj) + theme_minimal()

  ggsave(filename = paste0("plot/", dataset_name, "/mayo_matrix_", m, ".pdf"), dpi=400, plot = res_heatmap, width = 8, height = 8)
  ggsave(filename = paste0("plot/", dataset_name, "/volcano_", m, ".pdf"), dpi=400, plot = volc, width = 10, height = 8)
}

# test scMayoMap ####
comparison_tibble <- dplyr::tibble()
lfc_cut <- 1
for (m in c("devil", "nebula", "glmGamPoi")) {
  de_res <- de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()
  if (sum(grepl("ENSG", de_res$name)) == nrow(de_res)) {
    suppressMessages(de_res$name <- mapIds(org.Hs.eg.db, keys=de_res$name,column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
  }

  for (n_markers in c(5, 10, 25, 50, 100, 200, 300, 400, 500)) {
    print(n_markers)
    for (pval_cut in c(.05, 1e-5, 1e-10, 1e-25, 1e-50)) {
      print(pval_cut)
    #for (m in c("devil")) {

      # de_res$name <- mapIds(org.Hs.eg.db, keys=de_res$name,column="SYMBOL", keytype="ENSEMBL", multiVals="first")

      cluster_values <- de_res$cluster %>% unique()
      remove_genes <- grepl("^ENS", de_res$name)
      de_res <- de_res[!remove_genes, ]

      top_genes <- de_res %>%
        dplyr::group_by(cluster) %>%
        dplyr::filter(lfc > lfc_cut, adj_pval <= pval_cut) %>%
        dplyr::arrange(-lfc) %>%
        dplyr::slice(1:n_markers)

      c <- 1
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

      rez$true_score = lapply(1:nrow(rez), function(i) {
        r <- rez[i,]
        if (r$true_cell_type == r$pred) {
          return(r$score)
          #return(1)
        } else {
          return(0)
        }
      }) %>% unlist()
      acc = sum(rez$true_score) / length(unique(seurat_obj$seurat_clusters))

      comparison_tibble <- dplyr::bind_rows(
        comparison_tibble,
        dplyr::tibble(model=m, acc=acc, n_markers=n_markers, pval_cut=pval_cut)
      )
    }
  }
}

comparison_tibble %>%
  #dplyr::filter(pval_cut > 10) %>%
  #dplyr::filter(pval_cut <= 1e-25) %>%
  ggplot(mapping = aes(x=n_markers, y=acc, col=model)) +
  geom_point() +
  geom_smooth() +
  theme_bw() +
  scale_x_continuous(transform = 'log10') +
  labs(x = "N markers", y = "Accuracy", col="Algorithm")

ggsave(paste0("plot/", dataset_name, "/assigment_score.pdf"), dpi=300, width = 8, height = 5)

# # Extract average results ####
# n_markers <- 10
# pval_cut <- .05
# lfc_cut <- .2
# whole_results <- dplyr::tibble()
# #for (pval_cut in c(.05, 1e-5, 1e-10, 1e-25, 1e-50)) {
# for (pval_cut in c(.05, 1e-20, 1e-50)) {
#   print(pval_cut)
#   for (n_markers in c(5, 50, 100)) {
#     #for (n_markers in c(5, 100)) {
#     print(n_markers)
#     for (m in c("devil", "nebula", "glmGamPoi")) {
#       #for (m in c("devil", "nebula")) {
#       de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()
#
#       input_scMayo <- suppressMessages(
#         prepScMayoInput(de_res_total, as.matrix(seurat_obj@assays$RNA$counts), seurat_obj$seurat_clusters,
#                         n_markers = n_markers, lfc_cut = lfc_cut, pval_cut = pval_cut, distinct_marker = FALSE)
#       )
#
#       scMayoObj <- scMayoMap(data = input_scMayo, tissue = tissue, pct.cutoff = .25)
#       #seurat_obj$cell_type <- cell_type_names_to_scMayo_names(seurat_obj$cell_type, tissue)
#       ground_truth <- computeGroundTruth(seurat_obj)
#
#       pred_res <- lapply(1:nrow(scMayoObj$res), function(i) {
#         r <- scMayoObj$res[i,]
#         s <- strsplit(r$ES.norm, split = ";") %>% unlist() %>% as.numeric()
#         c <- strsplit(r$celltype, split = ";") %>% unlist()
#         c <- lapply(c, function(x) {stringr::str_replace_all(x, "cell", "")}) %>% unlist()
#         c <- lapply(c, function(x) {stringr::str_replace_all(x, " ", "")}) %>% unlist()
#         dplyr::tibble(cell_type_pred = c, score = s, cluster = r$cluster)
#       }) %>% do.call('bind_rows', .)
#
#       pred_res <- pred_res %>%
#         dplyr::group_by(cluster) %>%
#         dplyr::filter(score == max(score)) %>%
#         dplyr::mutate(n = n()) %>%
#         dplyr::mutate(cell_type_pred = ifelse(n == 1, cell_type_pred, "None")) %>%
#         dplyr::mutate(score = ifelse(n == 1, score, 0)) %>%
#         dplyr::distinct(cell_type_pred, score, cluster) %>%
#         dplyr::left_join(ground_truth, by='cluster')
#
#       pred_res$final_score <- lapply(1:nrow(pred_res), function(i) {
#         r <- pred_res[i,]
#         if (grepl(r$cell_type_pred, r$true_cell_type)) {
#           return(r$score)
#         } else {
#           return(0)
#         }
#       }) %>% unlist()
#
#       pred_res <- pred_res %>%
#         dplyr::select(cell_type_pred, cluster, true_cell_type, final_score) %>%
#         dplyr::mutate(model = m, pval_cut=pval_cut, n_markers=n_markers)
#
#       whole_results <- dplyr::bind_rows(whole_results, pred_res)
#     }
#   }
# }
#
# saveRDS(whole_results, paste0("results/", dataset_name, "_res_class.rds"))
#
# n_clusters <- length(unique(seurat_obj$seurat_clusters))
#
# whole_results %>%
#   #dplyr::filter(pval_cut <= 1e-50) %>%
#   group_by(model, pval_cut, n_markers) %>%
#   summarise(ACC = sum(final_score) / n_clusters) %>%
#   ggplot(mapping = aes(x=as.factor(n_markers), y=ACC, fill=model)) +
#   geom_violin() +
#   geom_jitter() +
#   theme_bw() +
#   labs(x = "N markers", y="Assignment score") +
#   theme(legend.position = 'bottom')
#
#
# whole_results %>%
#   #dplyr::filter(pval_cut <= 1e-50) %>%
#   group_by(model, pval_cut, n_markers) %>%
#   summarise(ACC = sum(final_score) / n_clusters) %>%
#   ggplot(mapping = aes(x=n_markers, y=ACC, fill=model, col=model)) +
#   geom_point() +
#   geom_smooth()
# geom_violin() +
#   geom_jitter() +
#   theme_bw() +
#   labs(x = "N markers", y="Assignment score") +
#   theme(legend.position = 'bottom')
#
# ggsave(paste0("plot/", dataset_name, "/assigment_score.pdf"), dpi=300, width = 8, height = 5)
#
# whole_results %>%
#   group_by(model, pval_cut, n_markers) %>%
#   summarise(IR = sum(final_score == 0) / n_clusters) %>%
#   ggplot(mapping = aes(x=as.factor(n_markers), y=IR, fill=model)) +
#   geom_violin() +
#   theme_bw() +
#   labs(x = "N markers", y="Non-assignment score") +
#   theme(legend.position = 'bottom')
# ggsave(paste0("plot/", dataset_name, "/indecisive_score.pdf"), dpi=300, width = 8, height = 5)
#
#
# de_devil <- readRDS(paste0('results/', dataset_name, '/', "devil", '.RDS')) %>% na.omit() %>% dplyr::filter(cluster == 4)
# de_nebula <- readRDS(paste0('results/', dataset_name, '/', "nebula", '.RDS')) %>% na.omit() %>% dplyr::filter(cluster == 4)
#
# gs = unique(c(de_devil$name, de_nebula$name))
# gs <- gs[(gs %in% de_devil$name) & (gs %in% de_nebula$name)]
# de_devil <- de_devil %>% dplyr::filter(name %in% gs)
# de_nebula <- de_nebula %>% dplyr::filter(name %in% gs)
#
# plot(-log10(de_devil$adj_pval), -log10(de_nebula$adj_pval))
# plot(de_devil$lfc, de_nebula$lfc)
#
#
#
#
