
rm(list = ls())
pkgs <- c("ggplot2", "dplyr","tidyr","tibble", "viridis", "smplot2", "Seurat", "VennDiagram", "gridExtra",
          "ggpubr", "ggrepel", "ggvenn", "ggpointdensity", "edgeR", "patchwork")
sapply(pkgs, require, character.only = TRUE)

cell_group_colors = c(
  "old" = "darkorange",
  "young" = "steelblue"
)

# input UMAPs ####
load("results/metadata_rna_umap.Rdata")

metadata_rna$cell_type %>% unique()

x_lims <- c(-5.5, 1)
y_lims <- c(0, 11)

box = dplyr::tibble(
  xmin=x_lims[1],
  xmax=x_lims[2],
  ymin=y_lims[1],
  ymax=y_lims[2],
)

umap_cell_types <- metadata_rna %>%
  #`rownames<-`(1:nrow(metadata_atac)) %>%
  dplyr::mutate(cell_type = if_else(cell_type %in% c("Myonuclei TI", "Myonuclei TII"), cell_type, "Other")) %>%
  #dplyr::filter(cell_type %in% c("Myonuclei TI", "Myonuclei TII")) %>%
  #dplyr::sample_n(5000) %>%
  dplyr::select(umap_1, umap_2, cell_type, age_pop) %>%
  dplyr::mutate(idx = row_number()) %>%
  #dplyr::filter(umap_1 >= -1, umap_2 >= -6, umap_2 <= 7.5) %>%
  ggplot(mapping = aes(x=umap_1, y=umap_2, col=cell_type)) +
  geom_point(alpha = .5, size = .5) +
  geom_rect(mapping = aes(xmin=x_lims[1], ymin=y_lims[1], xmax=x_lims[2], ymax=y_lims[2]), col="black", fill = alpha("white", 0)) +
  scale_color_manual(values = c("Myonuclei TI" = "maroon", "Myonuclei TII"="#880", "Other" = "gray90")) +
  theme_bw() +
  labs(x = "UMAP 1", y = "UMAP 2", col="Cell type")
umap_cell_types


metadata_rna %>%
  #`rownames<-`(1:nrow(metadata_atac)) %>%
  #dplyr::mutate(cell_type = if_else(cell_type %in% c("Myonuclei TI", "Myonuclei TII"), cell_type, "Other")) %>%
  dplyr::filter(cell_type %in% c("Myonuclei TI", "Myonuclei TII")) %>%
  #dplyr::sample_n(5000) %>%
  dplyr::select(umap_1, umap_2, cell_type, age_pop) %>%
  dplyr::mutate(idx = row_number()) %>%
  #dplyr::filter(umap_1 >= -1, umap_2 >= -6, umap_2 <= 7.5) %>%
  ggplot(mapping = aes(x=umap_1, y=umap_2, col=age_pop)) +
  geom_point(alpha = .5, size = .5) +
  theme_bw() +
  scale_color_manual(values = c("old_pop" = "goldenrod3", "young_pop" = "#483D8B")) +
  lims(x = x_lims) +
  labs(x = "UMAP 1", y = "UMAP 2", col="Cell type") +
  theme_void()









