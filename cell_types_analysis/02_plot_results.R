#setwd("~/GitHub/cell_types_analysis")
rm(list=ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble","reshape2", "AnnotationDbi", "org.Hs.eg.db")
sapply(pkgs, require, character.only = TRUE)
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

if (!(file.exists(paste0("results/", dataset_name)))) {
  dir.create(paste0("results/", dataset_name))
}

if (!(file.exists(paste0("plot/", dataset_name)))) {
  dir.create(paste0("plot/", dataset_name))
}

seurat_obj <- readRDS(paste0('results/', dataset_name, '_seurat.RDS'))

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
  x = seurat_obj@reductions$pca@cell.embeddings[,1],
  y = seurat_obj@reductions$pca@cell.embeddings[,2],
  donor_id = seurat_obj$donor,
  cell_type = seurat_obj$cell_type,
  cluster = seurat_obj$seurat_clusters,
)

saveRDS(d_umap, paste0("results/", dataset_name, "/umap_tibble.rds"))
ggsave(filename = paste0("plot/", dataset_name, "/umap_cluster.pdf"), dpi=400, plot = umap_plot_seurat, width = 8, height = 8)
ggsave(filename = paste0("plot/", dataset_name, "/umap_cell_types.pdf"), dpi=400, plot = umap_plot_labels, width = 8, height = 8)

m <- 'devil'
for (m in c("devil", "nebula", "glmGamPoi")) {
  de_res_total <- readRDS(paste0('results/', dataset_name, '_', m, '.RDS'))

  input_scMayo <- prepScMayoInput(
    de_res_total,
    as.matrix(seurat_obj@assays$RNA$counts),
    seurat_obj$seurat_clusters,
    n_markers = 50,
    lfc_cut = 1,
    pval_cut = .05,
    distinct_marker = FALSE) %>%
    dplyr::mutate(gene_id = paste(gene, cluster, sep = ":"))

  volc <- de_res_total %>%
    dplyr::mutate(name_id = paste(name, cluster, sep = ":")) %>%
    dplyr::mutate(is_marker = name_id %in% input_scMayo$gene_id) %>%
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

  obj <- scMayoMap(data = input_scMayo, tissue = tissue, pct.cutoff = 0)
  saveRDS(obj, paste0("results/", dataset_name, "/", m, "_scMayo.rds"))

  mayoMatrix <- plot_scMayoOutput(obj) + theme_minimal()

  ggsave(filename = paste0("plot/", dataset_name, "/mayo_matrix_", m, ".pdf"), dpi=400, plot = mayoMatrix, width = 8, height = 8)
  ggsave(filename = paste0("plot/", dataset_name, "/volcano_", m, ".pdf"), dpi=400, plot = volc, width = 10, height = 8)
}

# Extract average results ####
whole_results <- dplyr::tibble()
for (pval_cut in c(.05, .01, 1e-3, 1e-5,  1e-10,1e-20, 1e-40, 1e-50)) {
  print(pval_cut)
  for (n_markers in c(5, 10, 50, 100)) {
    print(n_markers)
    for (m in c("devil", "nebula", "glmGamPoi")) {
      de_res_total <- readRDS(paste0('results/', dataset_name, '_', m, '.RDS'))

      input_scMayo <- prepScMayoInput(de_res_total, as.matrix(seurat_obj@assays$RNA$counts), seurat_obj$seurat_clusters,
                                      n_markers = n_markers, lfc_cut = 1, pval_cut = pval_cut, distinct_marker = FALSE)

      scMayoObj <- scMayoMap(data = input_scMayo, tissue = tissue, pct.cutoff = 0)

      seurat_obj$cell_type <- cell_type_names_to_scMayo_names(seurat_obj$cell_type, tissue)
      ground_truth <- computeGroundTruth(seurat_obj)

      pred_res <- lapply(1:nrow(scMayoObj$res), function(i) {
        r <- scMayoObj$res[i,]
        s <- strsplit(r$ES.norm, split = ";") %>% unlist() %>% as.numeric()
        c <- strsplit(r$celltype, split = ";") %>% unlist()
        c <- lapply(c, function(x) {stringr::str_replace_all(x, "cell", "")}) %>% unlist()
        c <- lapply(c, function(x) {stringr::str_replace_all(x, " ", "")}) %>% unlist()
        dplyr::tibble(cell_type_pred = c, score = s, cluster = r$cluster)
      }) %>% do.call('bind_rows', .)

      pred_res <- pred_res %>%
        dplyr::group_by(cluster) %>%
        dplyr::filter(score == max(score)) %>%
        dplyr::mutate(n = n()) %>%
        dplyr::mutate(cell_type_pred = ifelse(n == 1, cell_type_pred, "None")) %>%
        dplyr::mutate(score = ifelse(n == 1, score, 0)) %>%
        dplyr::distinct(cell_type_pred, score, cluster) %>%
        dplyr::left_join(ground_truth, by='cluster')

      pred_res

      pred_res$final_score <- lapply(1:nrow(pred_res), function(i) {
        r <- pred_res[i,]
        if (grepl(r$cell_type_pred, r$true_cell_type)) {
          return(r$score)
        } else {
          return(0)
        }
      }) %>% unlist()

      pred_res <- pred_res %>%
        dplyr::select(cell_type_pred, cluster, true_cell_type, final_score) %>%
        dplyr::mutate(model = m, pval_cut=pval_cut, n_markers=n_markers)



      whole_results <- dplyr::bind_rows(whole_results, pred_res)
    }
  }
}

saveRDS(whole_results, paste0("results/", dataset_name, "_res_class.rds"))

n_clusters <- length(unique(seurat_obj$seurat_clusters))

whole_results %>%
  group_by(model, pval_cut, n_markers) %>%
  summarise(ACC = sum(final_score) / n_clusters) %>%
  ggplot(mapping = aes(x=as.factor(n_markers), y=ACC, fill=model)) +
  geom_violin() +
  geom_jitter() +
  theme_bw() +
  labs(x = "N markers", y="Assignment score") +
  theme(legend.position = 'bottom')

ggsave(paste0("plot/", dataset_name, "/assigment_score.pdf"), dpi=300, width = 8, height = 5)

whole_results %>%
  group_by(model, pval_cut, n_markers) %>%
  summarise(IR = sum(final_score == 0) / n_clusters) %>%
  ggplot(mapping = aes(x=as.factor(n_markers), y=IR, fill=model)) +
  geom_violin() +
  theme_bw() +
  labs(x = "N markers", y="Non-assignment score") +
  theme(legend.position = 'bottom')
ggsave(paste0("plot/", dataset_name, "/indecisive_score.pdf"), dpi=300, width = 8, height = 5)
