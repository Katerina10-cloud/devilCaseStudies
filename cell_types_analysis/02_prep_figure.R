rm(list=ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble","reshape2", "stringr", "AnnotationDbi", "scMayoMap", "Seurat", "org.Hs.eg.db")
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
  n.ct <- scTypeMapper %>% dplyr::filter(from == ct) %>% dplyr::distinct() %>% dplyr::pull(to)
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

lfc_cut <- 1
comparison_tibble <- dplyr::tibble()
m <- 'devil'
n_markers <- 25
pval_cut <- .01
for (m in c("devil", "nebula", "glmGamPoi")) {
  de_res <- de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()
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
  #dplyr::filter(pval_cut >= 1e-10) %>%
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

# # Comparison tibble our score ####
# anno <- scMayoMap::scMayoMapDatabase
# anno <- lapply(colnames(anno[2:ncol(anno)]), function(ct) {
#   dplyr::tibble(Type=ct, Marker = anno$gene[anno[,ct] == 1])
# }) %>% do.call('bind_rows', .) %>%
#   dplyr::filter(grepl(tissue, Type)) %>%
#   dplyr::mutate(Type = str_replace_all(Type, paste0(tissue, ":"), "")) %>%
#   dplyr::mutate(Type = str_replace_all(Type, paste0(" cell"), ""))
#
# lfc_cut <- 1
# m <- 'nebula'
# comparison_tibble <- dplyr::tibble()
# for (m in c("devil", "nebula", "glmGamPoi")) {
#   de_res <- de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()
#   if (sum(grepl("ENSG", de_res$name)) == nrow(de_res)) {
#     suppressMessages(de_res$name <- mapIds(org.Hs.eg.db, keys=de_res$name,column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
#   }
#   for (n_markers in c(5,10,25,50,100,500,1000)) {
#     print(n_markers)
#     #for (pval_cut in c(.05, 1e-5, 1e-10, 1e-20, 1e-30, 1e-40, 1e-50)) {
#     for (pval_cut in c(.05, .01, .001, 1e-10, 1e-20, 1e-50)) {
#       print(pval_cut)
#
#       de_res <- de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()
#       if (sum(grepl("ENSG", de_res$name)) == nrow(de_res)) {
#         suppressMessages(de_res$name <- mapIds(org.Hs.eg.db, keys=de_res$name,column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
#       }
#
#       cluster_values <- de_res$cluster %>% unique()
#       remove_genes <- grepl("^ENS", de_res$name)
#       de_res <- de_res[!remove_genes, ]
#
#       top_genes <- de_res %>%
#         dplyr::group_by(cluster) %>%
#         dplyr::filter(lfc > lfc_cut, adj_pval <= pval_cut) %>%
#         dplyr::arrange(-lfc) %>%
#         dplyr::slice(1:n_markers)
#
#       c <- 1
#       rez <- lapply(unique(top_genes$cluster), function(c) {
#         gg <- top_genes$name[top_genes$cluster == c]
#         table_pred <- anno %>%
#           dplyr::filter(Marker %in% gg) %>%
#           dplyr::pull(Type) %>%
#           table()
#         if (length(table_pred) == 0) {
#           dplyr::tibble(pred = NULL, score = 0, model=m, cluster=c)
#         } else {
#           dplyr::tibble(pred = names(table_pred), score = as.vector(table_pred) / sum(table_pred), model=m, cluster=c)
#         }
#       }) %>% do.call('bind_rows', .) %>%
#         dplyr::left_join(computeGroundTruth(seurat_obj), by='cluster') %>%
#         na.omit() %>%
#         dplyr::mutate(ground_truth = paste0(true_cell_type, " (", cluster, ")"))
#
#       rez <- rez %>%
#         dplyr::group_by(cluster) %>%
#         dplyr::filter(score == max(score)) %>%
#         dplyr::mutate(true_score = ifelse(true_cell_type == pred, score, 0))
#
#       acc = sum(rez$true_score) / length(unique(seurat_obj$seurat_clusters))
#       acc_1 = sum(rez$true_score > 0) / length(unique(seurat_obj$seurat_clusters))
#
#       comparison_tibble <- dplyr::bind_rows(
#         comparison_tibble,
#         dplyr::tibble(model=m, acc=acc, acc_1 =acc_1, n_markers=n_markers, pval_cut=pval_cut)
#       )
#     }
#   }
# }
#
# comparison_tibble %>%
#   #dplyr::filter(n_markers <= 100) %>%
#   dplyr::filter(pval_cut > 1e-03) %>%
#   dplyr::group_by(model, n_markers) %>%
#   dplyr::summarise(mean_acc = mean(acc), sd_acc = sd(acc)) %>%
#   #dplyr::summarise(mean_acc = mean(acc_1), sd_acc = sd(acc_1)) %>%
#   ggplot(mapping = aes(x=n_markers, y=mean_acc, ymin=mean_acc-sd_acc, ymax=mean_acc+sd_acc, col=model)) +
#   geom_pointrange(position=position_dodge(width=0.1)) +
#   geom_line(position=position_dodge(width=0.1)) +
#   theme_bw() +
#   scale_color_manual(values=method_colors) +
#   scale_x_continuous(transform = 'log10') +
#   labs(x = "N markers", y = "Classification score", col="Algorithm", fill="Algorithm")


# Volcano plot ####
if (dataset_name == "pbmc") {
  c <- 11
} else {
  c <- 1
}
m <- 'devil'

n_marker <- 5
lfc_cut <- 1
pval_adj <- 1e-50

de_res <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()

if (sum(grepl("ENSG", de_res$name)) == nrow(de_res)) {
  suppressMessages(de_res$name <- mapIds(org.Hs.eg.db, keys=de_res$name,column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
}
is.na(de_res$name) %>% sum()

cluster_values <- de_res$cluster %>% unique()
remove_genes <- is.na(de_res$name)
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

top_genes %>% dplyr::filter(cluster == c)

de_res %>%
  dplyr::filter(cluster == c) %>%
  dplyr::mutate(name_id = paste(name, cluster, sep = ":")) %>%
  dplyr::mutate(is_marker = name_id %in% top_genes$gene_id) %>%
  #dplyr::mutate(adj_pval = ifelse(adj_pval <= 6.181186e-16, 6.181186e-16, adj_pval)) %>%
  #dplyr::filter(!((!is_marker) & (adj_pval == 6.181186e-16))) %>%
  ggplot(mapping = aes(x=lfc, y=-log10(adj_pval), col=is_marker, size=-log10(adj_pval))) +
  geom_point() +
  scale_color_manual(values = c(alpha("black", .05), alpha(my_large_palette[c], 1))) +
  xlim(c(-5, NA)) +
  geom_hline(yintercept = -log10(pval_cut), linetype = "dashed") +
  geom_vline(xintercept = c(lfc_cut, -lfc_cut), linetype = "dashed") +
  theme_bw() +
  ggplot2::labs(x = expression(Log[2] ~ FC), y = expression(-log[10] ~ Pvalue), col="", size = "Significativity") +
  theme(
    text=element_text(size=12),
    legend.position = 'none'
  )
#theme(legend.position = 'bottom')

ggsave(filename = paste0(img_folder, "volcano_",c,".pdf"), dpi=300, width = 82, height = 85, units = 'mm')
if (save_svg) { ggsave(filename = paste0(img_folder, "volcano_",c,".svg"), dpi=300, width = 82, height = 85, units = 'mm') }

# Heatmap ####
lfc_cut <- 1
n_markers <- 50
pval_cut <- .05
m <- 'devil'

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

c <- 1
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
  unique(str_replace_all(v, " cell", ""))
}) %>% unlist()

rez_max <- rez %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(score == max(score)) %>%
  dplyr::mutate(is_correct = true_cell_type == pred)

# res_heatmap <- ggplot() +
#   # geom_point(rez_max %>% dplyr::filter(score > .1, is_correct), mapping = aes(x = ground_truth, y=pred, size=score * 3), shape=20, col="yellowgreen") +
#   # geom_point(rez_max %>% dplyr::filter(score > .1, !is_correct), mapping = aes(x = ground_truth, y=pred, size=score * 3), shape=20, col="indianred") +
#   geom_tile(rez %>% dplyr::filter(score > .1), mapping = aes(x=ground_truth, y=pred, fill=score), col='darkslategray') +
#   geom_point(rez_max %>% dplyr::filter(score > .1, is_correct), mapping = aes(x = ground_truth, y=pred), size =3, shape=20, col="yellowgreen") +
#   geom_point(rez_max %>% dplyr::filter(score > .1, !is_correct), mapping = aes(x = ground_truth, y=pred), shape=20, size=3, col="indianred") +
#   theme_bw() +
#   scale_fill_gradient(low = "#deebf7", high = "#08519c") +
#   labs(x = "True cell types", y = "Predicted cell types", fill="Score") +
#   #scale_color_continuous(type = "viridis") +
#   theme(
#     axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#     text=element_text(size=10),
#     legend.position = 'none'
#     )

res_heatmap <- ggplot() +
  # geom_point(rez_max %>% dplyr::filter(score > .1, is_correct), mapping = aes(x = ground_truth, y=pred, size=score * 3), shape=20, col="yellowgreen") +
  # geom_point(rez_max %>% dplyr::filter(score > .1, !is_correct), mapping = aes(x = ground_truth, y=pred, size=score * 3), shape=20, col="indianred") +
  geom_point(rez %>% dplyr::filter(score > .1), mapping = aes(x=true_cell_type, y=pred, fill=score, col=score, size=score * 5)) +
  geom_point(rez_max %>% dplyr::filter(score > .1, is_correct), mapping = aes(x = true_cell_type, y=pred, size=score), shape=20, col="yellowgreen") +
  geom_point(rez_max %>% dplyr::filter(score > .1, !is_correct), mapping = aes(x = true_cell_type, y=pred, size = score), shape=20, col="indianred") +
  theme_bw() +
  scale_color_gradient(low = "#deebf7", high = "#08519c") +
  scale_fill_gradient(low = "#deebf7", high = "#08519c") +
  labs(x = "True cell types", y = "Predicted cell types", fill="Score") +
  #scale_color_continuous(type = "viridis") +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    text=element_text(size=10),
    legend.position = 'none'
  )

res_heatmap

ggsave(filename = paste0(img_folder, "heatmap.pdf"), dpi=300, width = 140, height = 105, units = 'mm')
if (save_svg) { ggsave(filename = paste0(img_folder, "heatmap.svg"), dpi=300, width = 140, height = 105, units = 'mm') }

# Upset Plot ####
#mutations <- read.csv( system.file("extdata", "mutations.csv", package = "UpSetR"), header=T, sep = ",")

n_markers <- 1000000
pval_cut <- .05
all_markers <- c()
c <- 1

for (m in c("devil", "nebula", "glmGamPoi")) {
  de_res <- de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()
  if (sum(grepl("ENSG", de_res$name)) == nrow(de_res)) {
    suppressMessages(de_res$name <- mapIds(org.Hs.eg.db, keys=de_res$name,column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
  }

  top_genes <- de_res %>%
    dplyr::group_by(cluster) %>%
    dplyr::filter(lfc > lfc_cut, adj_pval <= pval_cut) %>%
    dplyr::arrange(-lfc) %>%
    dplyr::filter(cluster == c) %>%
    dplyr::slice(1:n_markers) %>% pull(name) %>% unname()

  all_markers <- c(all_markers, top_genes) %>% unique() %>% na.omit()
}

mm <- lapply(c("devil", "nebula", "glmGamPoi"), function(m) {
  de_res <- de_res_total <- readRDS(paste0('results/', dataset_name, '/', m, '.RDS')) %>% na.omit()
  if (sum(grepl("ENSG", de_res$name)) == nrow(de_res)) {
    suppressMessages(de_res$name <- mapIds(org.Hs.eg.db, keys=de_res$name,column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
  }

  top_genes <- de_res %>%
    dplyr::group_by(cluster) %>%
    dplyr::filter(lfc > lfc_cut, adj_pval <= pval_cut) %>%
    dplyr::arrange(-lfc) %>%
    dplyr::filter(cluster == c) %>%
    dplyr::slice(1:n_markers) %>% pull(name) %>% unname()

  c((all_markers %in% top_genes) %>% as.numeric())
}) %>% do.call('rbind', .)


mm_tibble <- lapply(1:ncol(mm), function(i) {
  print(i)
  g <- all_markers[i]
  v <- mm[,i] %>% unlist()
  dplyr::tibble(gene = g, devil=v[1], nebula=v[2], glmGamPoi=v[3])
}) %>% do.call("bind_rows", .) %>% as.data.frame()

method_order <- mm_tibble[,2:4] %>% colSums() %>% sort(decreasing = T)
method_colors_arranged <- c()

for (m in names(method_order)) {
  method_colors_arranged <- c(method_colors_arranged, method_colors[which(names(method_colors) == m)])
}

upset_plot <- UpSetR::upset(data.frame(mm_tibble), sets.bar.color = method_colors_arranged,
                            order.by = "freq", text.scale = 1.8)
svg(filename = paste0(img_folder, 'upset_plot.svg'), width = 9, height = 5)
upset_plot %>% print()
dev.off()

# Timing plot ####
timing <- readRDS(paste0("results/",dataset_name,"/time.RDS")) %>%
  dplyr::mutate(time = as.numeric(delta_time, units = "secs")) %>%
  dplyr::arrange(-time)

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
final_plot <- (pA / pCB) + plot_annotation(tag_levels = "A")
ggsave(paste0("plot_figure/", dataset_name, ".png"), dpi=400, width = 10, height = 7, plot = final_plot)
