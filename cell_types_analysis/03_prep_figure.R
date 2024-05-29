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
computeGroundTruth(seurat_obj) %>% dplyr::arrange(cluster)
seurat_obj$cell_type <- cell_type_names_to_scMayo_names(seurat_obj$cell_type, tissue)
computeGroundTruth(seurat_obj) %>% dplyr::arrange(cluster)

umap_plot_seurat <- Seurat::DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = T,
  repel = T) +
  ggtitle("") +
  scale_color_manual(values = my_large_palette) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

umap_plot_labels <- Seurat::DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "cell_type",
  label = T,
  repel = T) +
  ggtitle("") +
  scale_color_manual(values = my_large_palette) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

d_umap <- dplyr::tibble(
  x = seurat_obj@reductions$umap@cell.embeddings[,1],
  y = seurat_obj@reductions$umap@cell.embeddings[,2],
  donor_id = seurat_obj$donor,
  cell_type = seurat_obj$cell_type,
  cluster = seurat_obj$seurat_clusters,
)

ggsave(filename = paste0("plot_figure/umap_cluster.svg"), dpi=400, plot = umap_plot_seurat, width = 8, height = 8)
ggsave(filename = paste0("plot_figure/umap_labels.svg"), dpi=400, plot = umap_plot_labels, width = 8, height = 8)

anno <- scMayoMap::scMayoMapDatabase
anno <- lapply(colnames(anno[2:ncol(anno)]), function(ct) {
  dplyr::tibble(Type=ct, Marker = anno$gene[anno[,ct] == 1])
}) %>% do.call('bind_rows', .)
anno <- anno %>%
  dplyr::filter(grepl(tissue, Type)) %>%
  dplyr::mutate(Type = str_replace_all(Type, paste0(tissue, ":"), "")) %>%
  dplyr::mutate(Type = str_replace_all(Type, paste0(" cell"), ""))

lfc_cut <- 1
comparison_tibble <- dplyr::tibble()
for (m in c("devil", "nebula", "glmGamPoi")) {
  de_res <- de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()
  if (sum(grepl("ENSG", de_res$name)) == nrow(de_res)) {
    suppressMessages(de_res$name <- mapIds(org.Hs.eg.db, keys=de_res$name,column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
  }
  for (n_markers in c(5,10,20,50,100,200,300,1000)) {
    print(n_markers)
    #for (pval_cut in c(.05, 1e-5, 1e-10, 1e-20, 1e-30, 1e-40, 1e-50)) {
    for (pval_cut in c(.05, .01, .001, 1e-5, 1e-10, 1e-20, 1e-50)) {
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

      # rez$true_score = lapply(1:nrow(rez), function(i) {
      #   r <- rez[i,]
      #   if (r$true_cell_type == r$pred) {
      #     #return(r$score)
      #     return(1)
      #   } else {
      #     return(0)
      #   }
      # }) %>% unlist()
      acc = sum(rez$true_score) / length(unique(seurat_obj$seurat_clusters))

      comparison_tibble <- dplyr::bind_rows(
        comparison_tibble,
        dplyr::tibble(model=m, acc=acc, n_markers=n_markers, pval_cut=pval_cut)
      )
    }
  }
}

comparison_tibble %>%
  dplyr::filter(pval_cut > 1e-10) %>%
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

ggsave(paste0("plot_figure/assigment_score.svg"), dpi=300, width = 8, height = 5)


# Volcano plot and heatmap
m <- 'devil'
c <- 11
lfc_cut <- 1
pval_adj <- 1e-50

de_res <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()

if (sum(grepl("ENSG", de_res$name)) == nrow(de_res)) {
  suppressMessages(de_res$name <- mapIds(org.Hs.eg.db, keys=de_res$name,column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
}

cluster_values <- de_res$cluster %>% unique()
remove_genes <- grepl("^ENS", de_res$name)
de_res <- de_res[!remove_genes, ]

lfc_cut <- 1
n_markers <- 5
pval_cut <- .05

top_genes <- de_res %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(lfc > lfc_cut, adj_pval <= pval_cut) %>%
  dplyr::arrange(-lfc) %>%
  # dplyr::filter(abs(lfc) > lfc_cut, adj_pval <= pval_cut) %>%
  # dplyr::arrange(-abs(lfc)) %>%
  dplyr::slice(1:n_markers) %>%
  dplyr::mutate(gene_id = paste(name, cluster, sep=':'))

de_res %>%
  dplyr::filter(cluster == c) %>%
  dplyr::mutate(name_id = paste(name, cluster, sep = ":")) %>%
  dplyr::mutate(is_marker = name_id %in% top_genes$gene_id) %>%
  dplyr::mutate(adj_pval = ifelse(adj_pval <= 6.181186e-16, 6.181186e-16, adj_pval)) %>%
  dplyr::filter(!((!is_marker) & (adj_pval == 6.181186e-16))) %>%
  ggplot(mapping = aes(x=lfc, y=-log10(adj_pval), col=is_marker, size=1/adj_pval)) +
  geom_point() +
  scale_color_manual(values = c(alpha("black", .1), alpha(my_large_palette[c], 1))) +
  xlim(c(-5, NA)) +
  geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  geom_vline(xintercept = c(lfc_cut, -lfc_cut), linetype = "dashed") +
  theme_bw() +
  ggplot2::labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col="", size = "Significativity")
#theme(legend.position = 'bottom')

ggsave(filename = "plot_figure/volcano_11.svg", dpi=300, width = 8, height = 7)

lfc_cut <- 1
n_markers <- 1000
pval_cut <- .05

cluster_values <- de_res$cluster %>% unique()
remove_genes <- grepl("^ENS", de_res$name)
de_res <- de_res[!remove_genes, ]

top_genes <- de_res %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(lfc > lfc_cut, adj_pval <= pval_cut) %>%
  dplyr::arrange(-lfc) %>%
  dplyr::slice(1:n_markers)

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

rez %>%
  group_by(cluster) %>%
  dplyr::filter(score == max(score)) %>%
  dplyr::mutate(true_score = ifelse(true_cell_type == pred, 1, 0))

rez$true_score = lapply(1:nrow(rez), function(i) {
  r <- rez[i,]
  if (r$true_cell_type == r$pred) {
    #return(r$score)
    return(1)
  } else {
    return(0)
  }
}) %>% unlist()

acc = sum(rez$true_score) / length(unique(seurat_obj$seurat_clusters))
