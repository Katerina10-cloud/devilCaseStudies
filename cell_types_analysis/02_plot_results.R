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
tissue <- 'pancreas'

if (!(file.exists(paste0("results/", dataset_name)))) {
  dir.create(paste0("results/", dataset_name))
}

if (!(file.exists(paste0("plot/", dataset_name)))) {
  dir.create(paste0("plot/", dataset_name))
}

seurat_obj <- readRDS(paste0('results/', dataset_name, '/seurat.RDS'))
seurat_obj$cell_type <- cell_type_names_to_scMayo_names(seurat_obj$cell_type, tissue)

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
for (m in c("devil", "nebula", "glmGamPoi")) {
  de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()

  input_scMayo <- prepScMayoInput(
    de_res_total,
    as.matrix(seurat_obj@assays$RNA$counts),
    seurat_obj$seurat_clusters,
    n_markers = 10,
    lfc_cut = .2,
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
n_markers <- 10
pval_cut <- .05
lfc_cut <- .2
whole_results <- dplyr::tibble()
#for (pval_cut in c(.05, 1e-5, 1e-10, 1e-25, 1e-50)) {
for (pval_cut in c(.05, 1e-20, 1e-50)) {
  print(pval_cut)
  for (n_markers in c(5, 50, 100)) {
  #for (n_markers in c(5, 100)) {
    print(n_markers)
    for (m in c("devil", "nebula", "glmGamPoi")) {
    #for (m in c("devil", "nebula")) {
      de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()

      input_scMayo <- suppressMessages(
        prepScMayoInput(de_res_total, as.matrix(seurat_obj@assays$RNA$counts), seurat_obj$seurat_clusters,
                        n_markers = n_markers, lfc_cut = lfc_cut, pval_cut = pval_cut, distinct_marker = FALSE)
      )

      scMayoObj <- scMayoMap(data = input_scMayo, tissue = tissue, pct.cutoff = .25)
      #seurat_obj$cell_type <- cell_type_names_to_scMayo_names(seurat_obj$cell_type, tissue)
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
  #dplyr::filter(pval_cut <= 1e-50) %>%
  group_by(model, pval_cut, n_markers) %>%
  summarise(ACC = sum(final_score) / n_clusters) %>%
  ggplot(mapping = aes(x=as.factor(n_markers), y=ACC, fill=model)) +
  geom_violin() +
  geom_jitter() +
  theme_bw() +
  labs(x = "N markers", y="Assignment score") +
  theme(legend.position = 'bottom')


whole_results %>%
  #dplyr::filter(pval_cut <= 1e-50) %>%
  group_by(model, pval_cut, n_markers) %>%
  summarise(ACC = sum(final_score) / n_clusters) %>%
  ggplot(mapping = aes(x=n_markers, y=ACC, fill=model, col=model)) +
  geom_point() +
  geom_smooth()
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



# test scSorter
library(scSorter)
tissue <- 'pancreas'
lfc_cut <- .2

scSorter_res <- dplyr::tibble()
for (m in c('devil', 'nebula', 'glmGamPoi')) {
  print(m)
  for (pval_cut in c(.05, .01, 1e-10, 1e-25, 1e-50)) {
    print(paste0("  ", pval_cut))
    for (n_markers in c(5, 10, 25, 50, 100)) {
      print(paste0("    ", n_markers))
      de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()
      gs <- lapply(unique(seurat_obj$seurat_clusters), function(c) {
        de_res_total %>%
          dplyr::filter(cluster == c, lfc > lfc_cut, adj_pval <= pval_cut) %>%
          dplyr::arrange(-lfc) %>%
          dplyr::slice(1:n_markers)
      }) %>% do.call("bind_rows", .) %>% pull(name) %>% unique()

      expr <- seurat_obj@assays$RNA$counts

      expr <- lapply(unique(seurat_obj$seurat_clusters), function(c) {
        expr[,seurat_obj$seurat_clusters == c] %>% rowSums()
      }) %>% do.call("cbind", .)

      expr = xnormalize_scData(expr)
      expr <- expr[gs, ]

      anno <- scMayoMap::scMayoMapDatabase
      anno <- lapply(colnames(anno[2:ncol(anno)]), function(ct) {
        dplyr::tibble(Type=ct, Marker = anno$gene[anno[,ct] == 1])
      }) %>% do.call('bind_rows', .)
      anno <- anno %>%
        dplyr::filter(grepl(tissue, Type)) %>%
        # dplyr::group_by(Marker) %>%
        # dplyr::mutate(n = n()) %>%
        # dplyr::filter(n == 1) %>%
        # dplyr::select(Type, Marker) %>%
        dplyr::filter(Marker %in% gs)

      cell_classification <- scSorter::scSorter(expr, as.data.frame(anno), default_weight = 2)

      cell_classification$Pred_Type
      avg_score <- computeGroundTruth(seurat_obj) %>%
        dplyr::mutate(pred = cell_classification$Pred_Type) %>%
        dplyr::group_by(cluster) %>%
        dplyr::mutate(score = grepl(tolower(true_cell_type), tolower(pred))) %>%
        pull(score) %>% mean()

      scSorter_res <- dplyr::bind_rows(scSorter_res, dplyr::tibble(pval_cut=pval_cut, n_markers=n_markers, model=m, score=avg_score))
    }
  }
}

scSorter_res %>%
  ggplot(mapping = aes(x = n_markers, y=score, col=model, fill=model)) +
  geom_point() +
  geom_smooth() +
  scale_x_continuous(transform = "log10")

scSorter_res %>%
  #dplyr::filter(pval_cut == 1e-20) %>%
  ggplot(mapping = aes(x = as.factor(n_markers), y=score, col=model)) +
  geom_violin()
  #scale_x_continuous(transform = "log10")


de_devil <- readRDS(paste0('results/', dataset_name, '/', "devil", '.RDS')) %>% na.omit() %>% dplyr::filter(cluster == 4)
de_nebula <- readRDS(paste0('results/', dataset_name, '/', "nebula", '.RDS')) %>% na.omit() %>% dplyr::filter(cluster == 4)

gs = unique(c(de_devil$name, de_nebula$name))
gs <- gs[(gs %in% de_devil$name) & (gs %in% de_nebula$name)]
de_devil <- de_devil %>% dplyr::filter(name %in% gs)
de_nebula <- de_nebula %>% dplyr::filter(name %in% gs)

plot(-log10(de_devil$adj_pval), -log10(de_nebula$adj_pval))
plot(de_devil$lfc, de_nebula$lfc)

