rm(list=ls())
pkgs <- c("patchwork", "ggpubr","ggplot2", "dplyr","tidyr","tibble","reshape2", "stringr", "AnnotationDbi", "scMayoMap", "Seurat", "org.Hs.eg.db")
sapply(pkgs, require, character.only = TRUE)
source("utils.R")

set.seed(SEED)

args = commandArgs(trailingOnly=TRUE)

## Input data
dataset_name <- args[1]
tissue <- args[2]
save_svg <- as.logical(args[3])

dataset_name <- "BaronPancreasData"
tissue <- "pancreas"
save_svg <- F

img_folder <- paste0("plot_figure/", dataset_name, "/")
if (!dir.exists(img_folder)) {
  dir.create(img_folder)
}

if (!(file.exists(paste0("results/", dataset_name)))) {
  dir.create(paste0("results/", dataset_name))
}

scTypeMapper <- read.delim("scTypeMapper.csv", sep = ",") %>% dplyr::rename(Tissue = tissue) %>% dplyr::filter(Tissue == tissue)

seurat_obj <- readRDS(paste0('results/', dataset_name, '/seurat.RDS'))
computeGroundTruth(seurat_obj) %>% dplyr::arrange(cluster)
print("Seurat obj cell types")
print(seurat_obj$cell_type %>% unique())


anno <- scMayoMap::scMayoMapDatabase
anno <- lapply(colnames(anno[2:ncol(anno)]), function(ct) {
  dplyr::tibble(Type=ct, Marker = anno$gene[anno[,ct] == 1])
}) %>% do.call('bind_rows', .)
anno <- anno %>%
  dplyr::filter(grepl(tissue, Type)) %>%
  dplyr::mutate(Type = str_replace_all(Type, paste0(tissue, ":"), "")) %>%
  dplyr::mutate(Type = str_replace_all(Type, paste0(" cell"), ""))
print("ScMayo cell types")
print(anno$Type %>% unique())

new_cell_type <- rep(NA, length(seurat_obj$cell_type))
for (ct in unique(seurat_obj$cell_type)) {
  print(ct)
  idx <- which(seurat_obj$cell_type == ct)
  n.ct <- scTypeMapper %>% dplyr::filter(from == ct) %>% dplyr::distinct() %>% dplyr::pull(to) %>% unique()
  new_cell_type[idx] <- n.ct
}
seurat_obj$cell_type <- new_cell_type
ground_truth <- computeGroundTruth(seurat_obj) %>% dplyr::arrange(cluster)
# seurat_obj$cell_type <- cell_type_names_to_scMayo_names(seurat_obj$cell_type, tissue)
# computeGroundTruth(seurat_obj) %>% dplyr::arrange(cluster)

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

data_umap <- dplyr::tibble(
  umap1=seurat_obj@reductions$umap@cell.embeddings[,1],
  umap2=seurat_obj@reductions$umap@cell.embeddings[,2],
  cluster=seurat_obj$seurat_clusters,
  cell_type=seurat_obj$cell_type
) %>%
  dplyr::left_join(ground_truth, by = "cluster") %>%
  dplyr::mutate(legend = paste(cluster, ". ", true_cell_type))

strings <- unique(data_umap$legend)
extract_numbers <- function(x) {
  as.numeric(sub("\\..*", "", x))
}
numbers <- sapply(strings, extract_numbers)
sorted_strings <- strings[order(numbers)]
data_umap$legend = factor(data_umap$legend, levels = sorted_strings)

data_avg_umap <- data_umap %>%
  dplyr::group_by(cluster) %>%
  dplyr::select(umap1, umap2) %>%
  dplyr::summarise_all(mean)

pA <- ggplot() +
  geom_point(data_umap, mapping = aes(x=umap1, y=umap2, col=legend), size=.1) +
  geom_point(data = data_avg_umap, mapping = aes(x=umap1, y = umap2)) +
  scale_color_manual(values = my_large_palette) +
  theme_bw() +
  labs(x = "UMAP 1", y = "UMAP 2", col="Cluster") +
  guides(color = guide_legend(override.aes = list(size=2))) +
  geom_label(data_avg_umap, mapping=aes(x=umap1, y = umap2, label = cluster)) +
  theme(text=element_text(size=12))
pA
# umap_plot_labels <- Seurat::DimPlot(
#   seurat_obj,
#   reduction = "umap",
#   group.by = "cell_type",
#   label = T,
#   repel = T) +
#   ggtitle("") +
#   scale_color_manual(values = my_large_palette) +
#   theme(axis.line=element_blank(),axis.text.x=element_blank(),
#         axis.text.y=element_blank(),axis.ticks=element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),legend.position="none",
#         panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),plot.background=element_blank())

ggsave(filename = paste0(img_folder, "umap_cluster.pdf"), dpi=400, plot = umap_plot_seurat, width = 8, height = 8)
if (save_svg) { ggsave(filename = paste0(img_folder, "umap_cluster.svg"), dpi=400, plot = umap_plot_seurat, width = 8, height = 8) }
#ggsave(filename = paste0("plot_figure/umap_labels.svg"), dpi=400, plot = umap_plot_labels, width = 8, height = 8)

# anno <- scMayoMap::scMayoMapDatabase
# anno <- lapply(colnames(anno[2:ncol(anno)]), function(ct) {
#   dplyr::tibble(Type=ct, Marker = anno$gene[anno[,ct] == 1])
# }) %>% do.call('bind_rows', .)
# anno <- anno %>%
#   dplyr::filter(grepl(tissue, Type)) %>%
#   dplyr::mutate(Type = str_replace_all(Type, paste0(tissue, ":"), "")) %>%
#   dplyr::mutate(Type = str_replace_all(Type, paste0(" cell"), ""))
# anno$Type %>% unique()
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
rm(counts)

# Comparison tibble ####
#lfc_cut <- 1
#pval_cut <- .05
comparison_tibble <- dplyr::tibble()
for (m in c("devil", "nebula", "glmGamPoi")) {
  de_res <- de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()
  if (m == "nebula") {
    de_res$lfc <- de_res$lfc / log(2)
  }
  if (sum(grepl("ENSG", de_res$name)) == nrow(de_res)) {
    suppressMessages(de_res$name <- mapIds(org.Hs.eg.db, keys=de_res$name,column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
  }
  for (n_markers in c(5,10,25,50,100,500,1000)) {
    print(n_markers)
    #for (pval_cut in c(.05, 1e-5, 1e-10, 1e-20, 1e-30, 1e-40, 1e-50)) {
    for (pval_cut in c(.05, .01, .001, 1e-20)) {
      print(pval_cut)

      de_res <- de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()
      if (sum(grepl("ENSG", de_res$name)) == nrow(de_res)) {
        suppressMessages(de_res$name <- mapIds(org.Hs.eg.db, keys=de_res$name,column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
      }

      cluster_values <- de_res$cluster %>% unique()
      remove_genes <- is.na(de_res$name)
      de_res <- de_res[!remove_genes, ]

      top_genes <- de_res %>%
        dplyr::group_by(cluster) %>%
        dplyr::filter(lfc > lfc_cut, adj_pval <= pval_cut) %>%
        dplyr::arrange(-lfc) %>%
        dplyr::slice(1:n_markers)

      if (nrow(top_genes) > 0) {

        scMayoInput <- lapply(unique(top_genes$cluster), function(c) {
          top_genes %>%
            dplyr::filter(cluster == c) %>% pull(name) %>% table() %>% max()

          top_genes %>%
            dplyr::filter(cluster == c) %>%
            dplyr::left_join(percentage_tibble %>% dplyr::filter(cluster == c) %>% dplyr::select(!cluster), by='name') %>%
            dplyr::ungroup() %>%
            dplyr::select(name, pval, adj_pval, lfc, pct.1, pct.2, cluster) %>%
            dplyr::rename(p_val_adj = adj_pval, gene=name, avg_log2FC=lfc)
        }) %>% do.call('bind_rows', .)

        scMayoRes <- scMayoMap::scMayoMap(scMayoInput, scMayoMap::scMayoMapDatabase, pct.cutoff = 0, tissue = tissue)
        rez <- scMayoRes$markers %>%
          dplyr::left_join(computeGroundTruth(seurat_obj), by='cluster') %>%
          dplyr::mutate(pred = str_replace_all(celltype, " cell", "")) %>%
          dplyr::mutate(score = as.numeric(score)) %>%
          dplyr::filter(score > 0)
        rez$pred <- lapply(rez$pred, function(ct) {
          v <- scTypeMapper %>% dplyr::filter(from == ct) %>% dplyr::distinct() %>%pull(to)
          unique(v)
        }) %>% unlist()

        rez <- rez %>%
          dplyr::select(cluster, score, true_cell_type, pred) %>%
          dplyr::group_by(cluster, pred) %>%
          dplyr::mutate(score = sum(score)) %>%
          dplyr::ungroup() %>%
          dplyr::distinct() %>%
          dplyr::group_by(cluster) %>%
          dplyr::filter(score == max(score)) %>%
          dplyr::mutate(true_score = ifelse(true_cell_type == pred, score, 0))

        n_max <- length(unique(seurat_obj$seurat_clusters))

        acc = sum(rez$true_score) / n_max
        acc_1 = sum(rez$true_score > 0) / n_max

        comparison_tibble <- dplyr::bind_rows(
          comparison_tibble,
          dplyr::tibble(model=m, acc=acc, acc_1 =acc_1, n_markers=n_markers, pval_cut=pval_cut)
        )
      }
    }
  }
}

pB <- comparison_tibble %>%
  #dplyr::filter(n_markers <= 100) %>%
  dplyr::filter(pval_cut >= 1e-10) %>%
  dplyr::group_by(model, n_markers) %>%
  dplyr::summarise(mean_acc = mean(acc), sd_acc = sd(acc)) %>%
  #dplyr::summarise(mean_acc = mean(acc_1), sd_acc = sd(acc_1)) %>%
  ggplot(mapping = aes(x=n_markers, y=mean_acc, ymin=mean_acc-sd_acc, ymax=mean_acc+sd_acc, col=model)) +
  geom_pointrange(position=position_dodge(width=0.1)) +
  geom_line(position=position_dodge(width=0.1)) +
  theme_bw() +
  scale_color_manual(values=method_colors) +
  scale_x_continuous(trans = 'log10') +
  labs(x = "N markers", y = "Classification score", col="Algorithm", fill="Algorithm") +
  theme(
    text=element_text(size=12),
    legend.position = 'right'
  )
pB

ggsave(filename = paste0(img_folder, "assigment_score.pdf"), plot = last_plot(), dpi=300, width = 150, height = 60, units = "mm")
if (save_svg) { ggsave(filename = paste0(img_folder, "assigment_score.svg"), plot = last_plot(), dpi=300, width = 150, height = 60, units = "mm") }

# Gene overlap ####
lfc_cut <- 1
pval_cut <- .05
m = "nebula"
de_genes_called <- lapply(c("devil", "nebula", "glmGamPoi"), function(m) {
  de_res <- de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()
  if (m == "nebula") {
    de_res$lfc <- de_res$lfc / log(2)
  }
  de_res %>%
    dplyr::filter(lfc >= lfc_cut, adj_pval <= pval_cut) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::distinct(cluster, n) %>%
    dplyr::mutate(method = m)
}) %>% do.call('bind_rows', .)

de_genes_called %>%
  ggplot(mapping = aes(x = as.numeric(cluster), y=n, col=method)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = method_colors) +
  theme_bw()

# Timing plot ####
timing <- readRDS(paste0("results/",dataset_name,"/time.RDS")) %>%
  dplyr::mutate(time = as.numeric(delta_time, units = "secs")) %>%
  dplyr::arrange(-time)

timing %>%
  dplyr::left_join(de_genes_called, by="method") %>%
  dplyr::mutate(time_ratio = time / time[method == "devil"]) %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(genes_ratio = n / n[method == "devil"]) %>%
  dplyr::filter(method != "devil") %>%
  ggplot(mapping = aes(x = time_ratio, y = genes_ratio, col = method)) +
  geom_point()


pC <- timing %>%
  ggplot(mapping = aes(x=method, y=time, fill=method)) +
  geom_col() +
  theme_bw() +
  labs(x = "", y="Time (s)") +
  scale_fill_manual(values = method_colors) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    text=element_text(size=12),
    legend.position = 'none'
  )
pC

ggsave(filename = paste0(img_folder, "timing.pdf"), dpi=300, width = 40, height = 60, units = 'mm')
if (save_svg) {ggsave(filename = paste0(img_folder, "timing.svg"), dpi=300, width = 40, height = 60, units = 'mm')}

unlink("Rplots.pdf")

# save figure

pCB <- (pC | pB) + plot_layout(widths = c(1.2,3))
final_plot <- (pA / pCB) + plot_annotation(tag_levels = "A") + plot_layout(heights = c(2,1))
ggsave(paste0("plot_figure/", dataset_name, ".png"), dpi=400, width = 8, height = 8, plot = final_plot)
